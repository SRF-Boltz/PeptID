import datetime
import hashlib
import os
import pickle
import re
import tempfile
import time
import urllib.error
import urllib.request
import zipfile
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed

import pandas as pd


# ============================================================
# PeptID paths and downloadable library configuration
# ============================================================

# PeptID uses the directory containing this script as its root directory.
base_path = os.path.dirname(os.path.abspath(__file__))

results_path = os.path.join(base_path, "Results")
libraries_path = os.path.join(base_path, "Libraries")
queries_path = os.path.join(base_path, "Queries")

os.makedirs(results_path, exist_ok=True)
os.makedirs(libraries_path, exist_ok=True)
os.makedirs(queries_path, exist_ok=True)

GITHUB_USER = "SRF-Boltz"
GITHUB_REPO = "PeptID"

# These filenames must exactly match the assets attached to the latest
# published GitHub release.
LIBRARY_ZIPS = (
    "Function_Library.zip",
    "Source_Library.zip",
    "Combined_Library.zip",
)

COMBINED_LIBRARY_FILENAME = "Combined_Library.txt"
LIBRARY_INSTALL_MARKER = ".peptid_libraries_installed"

# Leave as None to have PeptID automatically use the single FASTA-like file
# placed in Queries. Alternatively, set this to a specific filename.
QUERY_FILENAME = None
QUERY_EXTENSIONS = (".fasta", ".fa", ".faa", ".fas", ".fna", ".txt")

MAX_WORKERS = 1
CHUNK_SIZE = 10
MIN_PEPTIDE_LENGTH = 5
CASE_SENSITIVE = False
VERBOSE_EVERY = 60
SAVE_EVERY_COMPLETED = 25

# If True, ignores all peptide entries containing X/x wildcard residues.
IGNORE_WILDCARD_PEPTIDES = True

# If True, excludes peptides where only one position is a real amino-acid
# residue and all remaining positions are X/x wildcards.
EXCLUDE_SINGLE_RESIDUE_WILDCARD_PEPTIDES = True

# If True, shorter peptide matches fully contained within a longer match
# in the same query sequence are removed from the final output.
KEEP_ONLY_LONGEST_OVERLAPPING_MATCHES = True


# ============================================================
# First-run setup: download and extract GitHub release libraries
# ============================================================

def _download_progress(block_count, block_size, total_size):
    downloaded = block_count * block_size
    if total_size > 0:
        downloaded = min(downloaded, total_size)
        percent = downloaded * 100 / total_size
        print(
            f"\r  {downloaded / (1024 ** 2):.1f} / "
            f"{total_size / (1024 ** 2):.1f} MiB ({percent:.1f}%)",
            end="",
            flush=True,
        )
    else:
        print(
            f"\r  {downloaded / (1024 ** 2):.1f} MiB downloaded",
            end="",
            flush=True,
        )


def _zip_is_valid(zip_path):
    if not os.path.isfile(zip_path):
        return False

    try:
        if not zipfile.is_zipfile(zip_path):
            return False
        with zipfile.ZipFile(zip_path, "r") as archive:
            archive.infolist()
        return True
    except (OSError, zipfile.BadZipFile):
        return False


def _download_release_asset(zip_name, zip_path):
    url = (
        f"https://github.com/{GITHUB_USER}/{GITHUB_REPO}/"
        f"releases/latest/download/{zip_name}"
    )
    partial_path = f"{zip_path}.part"

    if os.path.exists(partial_path):
        os.remove(partial_path)

    print(f"Downloading {zip_name} from the latest PeptID release...")

    try:
        urllib.request.urlretrieve(
            url,
            partial_path,
            reporthook=_download_progress,
        )
        print()
        os.replace(partial_path, zip_path)
    except urllib.error.HTTPError as exc:
        if os.path.exists(partial_path):
            os.remove(partial_path)
        raise RuntimeError(
            f"Could not download {zip_name} (HTTP {exc.code}). "
            "Confirm that the SRF-Boltz/PeptID repository is public, a release "
            "has been published, and the release asset name matches exactly."
        ) from exc
    except urllib.error.URLError as exc:
        if os.path.exists(partial_path):
            os.remove(partial_path)
        raise RuntimeError(
            f"Could not download {zip_name}. Check the internet connection. "
            f"Details: {exc.reason}"
        ) from exc
    except OSError as exc:
        if os.path.exists(partial_path):
            os.remove(partial_path)
        raise RuntimeError(
            f"Could not save {zip_name} in {libraries_path}: {exc}"
        ) from exc

    if not _zip_is_valid(zip_path):
        raise RuntimeError(
            f"The downloaded file is not a valid ZIP archive: {zip_path}"
        )


def _safe_extract_zip(zip_path, destination):
    destination_real = os.path.realpath(destination)

    with zipfile.ZipFile(zip_path, "r") as archive:
        for member in archive.infolist():
            member_path = os.path.realpath(
                os.path.join(destination_real, member.filename)
            )
            try:
                inside_destination = (
                    os.path.commonpath([destination_real, member_path])
                    == destination_real
                )
            except ValueError:
                inside_destination = False

            if not inside_destination:
                raise RuntimeError(
                    f"Unsafe path found in {os.path.basename(zip_path)}: "
                    f"{member.filename}"
                )

        archive.extractall(destination_real)


def _find_combined_library():
    direct_path = os.path.join(libraries_path, COMBINED_LIBRARY_FILENAME)
    if os.path.isfile(direct_path):
        return direct_path

    matches = []
    for root, _, files in os.walk(libraries_path):
        if COMBINED_LIBRARY_FILENAME in files:
            matches.append(os.path.join(root, COMBINED_LIBRARY_FILENAME))

    if len(matches) == 1:
        return matches[0]
    if len(matches) > 1:
        raise RuntimeError(
            "More than one Combined_Library.txt was found after extraction: "
            + ", ".join(matches)
        )
    return None


def install_libraries_if_needed():
    marker_path = os.path.join(libraries_path, LIBRARY_INSTALL_MARKER)
    combined_library = _find_combined_library()
    zip_paths = [os.path.join(libraries_path, name) for name in LIBRARY_ZIPS]

    installation_complete = (
        combined_library is not None
        and os.path.isfile(marker_path)
        and all(_zip_is_valid(path) for path in zip_paths)
    )

    if installation_complete:
        print("PeptID libraries found; download skipped.")
        return combined_library

    print("Installing PeptID libraries...")

    for zip_name, zip_path in zip(LIBRARY_ZIPS, zip_paths):
        if _zip_is_valid(zip_path):
            print(f"Using existing archive: {zip_name}")
        else:
            if os.path.exists(zip_path):
                os.remove(zip_path)
            _download_release_asset(zip_name, zip_path)

        print(f"Extracting {zip_name}...")
        _safe_extract_zip(zip_path, libraries_path)

    combined_library = _find_combined_library()
    if combined_library is None:
        raise FileNotFoundError(
            f"{COMBINED_LIBRARY_FILENAME} was not found after extracting "
            "Function_Library.zip, Source_Library.zip, and "
            "Combined_Library.zip. Check the contents of the release assets."
        )

    with open(marker_path, "w", encoding="utf-8") as marker:
        marker.write(
            "PeptID libraries installed from the latest GitHub release.\n"
            f"Repository: {GITHUB_USER}/{GITHUB_REPO}\n"
            f"Installed: {datetime.datetime.now().isoformat(timespec='seconds')}\n"
        )

    print("PeptID library installation complete.")
    return combined_library


def select_query_file():
    if QUERY_FILENAME:
        requested_path = os.path.join(queries_path, QUERY_FILENAME)
        if not os.path.isfile(requested_path):
            raise FileNotFoundError(
                f"Configured query FASTA file not found: {requested_path}"
            )
        return requested_path

    candidates = sorted(
        os.path.join(queries_path, filename)
        for filename in os.listdir(queries_path)
        if os.path.isfile(os.path.join(queries_path, filename))
        and filename.lower().endswith(QUERY_EXTENSIONS)
        and not filename.startswith(".")
    )

    if not candidates:
        raise FileNotFoundError(
            "No query FASTA file was found. Place one FASTA-like file "
            f"({', '.join(QUERY_EXTENSIONS)}) in: {queries_path}"
        )
    if len(candidates) > 1:
        candidate_names = ", ".join(os.path.basename(path) for path in candidates)
        raise RuntimeError(
            "Multiple query files were found in the Queries folder. Leave only "
            f"the file to analyse, or set QUERY_FILENAME. Found: {candidate_names}"
        )

    return candidates[0]


# ============================================================
# Global worker state
# ============================================================

PEPTIDE_INDEX = None


def init_worker(peptide_index):
    global PEPTIDE_INDEX
    PEPTIDE_INDEX = peptide_index


# ============================================================
# FASTA parsing
# ============================================================

def parse_fasta_to_df(filepath):
    records = []
    current_id = None
    current_seq = []

    with open(filepath, "r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue

            if line.startswith(">"):
                if current_id is not None:
                    records.append((current_id, "".join(current_seq)))
                current_id = line[1:].strip()
                current_seq = []
            else:
                current_seq.append(line)

    if current_id is not None:
        records.append((current_id, "".join(current_seq)))

    return pd.DataFrame(records, columns=["Q1", "Q2"])


# ============================================================
# Hashing and progress I/O
# ============================================================

def compute_file_hash(filepath, algorithm="sha256"):
    hasher = hashlib.new(algorithm)
    with open(filepath, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            hasher.update(chunk)
    return hasher.hexdigest()


def load_progress(progress_path):
    if progress_path and os.path.exists(progress_path):
        with open(progress_path, "rb") as handle:
            return pickle.load(handle)
    return None


def save_progress_atomic(progress_path, progress_data):
    if not progress_path:
        return

    os.makedirs(os.path.dirname(progress_path), exist_ok=True)

    fd, tmp_path = tempfile.mkstemp(
        prefix=".progress_",
        suffix=".tmp",
        dir=os.path.dirname(progress_path)
    )

    try:
        with os.fdopen(fd, "wb") as handle:
            pickle.dump(progress_data, handle, protocol=pickle.HIGHEST_PROTOCOL)
        os.replace(tmp_path, progress_path)
    finally:
        if os.path.exists(tmp_path):
            os.remove(tmp_path)


# ============================================================
# Peptide preprocessing
# ============================================================

def peptide_to_regex(peptide):
    pattern = "".join("." if aa.upper() == "X" else re.escape(aa) for aa in peptide)
    return re.compile(f"^{pattern}$")


def count_non_wildcard_residues(peptide):
    return sum(1 for aa in peptide if aa.upper() != "X")


def build_peptide_index_from_file(
    library_file,
    case_sensitive=False,
    min_peptide_length=5,
    ignore_wildcard_peptides=False,
    exclude_single_residue_wildcard_peptides=True
):
    exact_by_length = defaultdict(lambda: defaultdict(list))
    degenerate_by_length_and_first = defaultdict(lambda: defaultdict(list))
    degenerate_by_length_any_first = defaultdict(list)

    total_entries = 0
    indexed_entries = 0
    malformed_removed = 0
    short_removed = 0
    wildcard_removed = 0
    single_wildcard_removed = 0

    with open(library_file, "r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            total_entries += 1
            line = line.rstrip("\n\r")

            parts = line.split("\t")

            if len(parts) < 3:
                malformed_removed += 1
                continue

            peptide_id = str(parts[0])
            peptide_source = str(parts[1])
            peptide_original = str(parts[2]).strip()

            if not peptide_original:
                malformed_removed += 1
                continue

            peptide_match = peptide_original if case_sensitive else peptide_original.upper()
            contains_wildcard = "X" in peptide_match.upper()

            if len(peptide_match) < min_peptide_length:
                short_removed += 1
                continue

            if ignore_wildcard_peptides and contains_wildcard:
                wildcard_removed += 1
                continue

            if contains_wildcard and exclude_single_residue_wildcard_peptides:
                if count_non_wildcard_residues(peptide_match) <= 1:
                    single_wildcard_removed += 1
                    continue

            length = len(peptide_match)

            metadata = {
                "peptide_id": peptide_id,
                "peptide_source": peptide_source,
                "peptide_sequence": peptide_original,
                "peptide_sequence_match": peptide_match,
                "length": length
            }

            if contains_wildcard:
                metadata["regex"] = peptide_to_regex(peptide_match)

                first_residue = peptide_match[0].upper()
                if first_residue == "X":
                    degenerate_by_length_any_first[length].append(metadata)
                else:
                    degenerate_by_length_and_first[length][first_residue].append(metadata)
            else:
                exact_by_length[length][peptide_match].append(metadata)

            indexed_entries += 1

    print(f"Loaded peptide library with {total_entries} peptide entries")
    print(f"Removed malformed/empty entries: {malformed_removed}")
    print(f"Removed peptides shorter than {min_peptide_length} residues: {short_removed}")

    if ignore_wildcard_peptides:
        print(f"Ignored wildcard-containing peptides: {wildcard_removed}")
    elif exclude_single_residue_wildcard_peptides:
        print(
            f"Excluded {single_wildcard_removed} peptides containing only one "
            f"non-wildcard residue and otherwise X/x wildcard positions"
        )

    return {
        "exact_by_length": dict(exact_by_length),
        "degenerate_by_length_and_first": {
            length: dict(first_groups)
            for length, first_groups in degenerate_by_length_and_first.items()
        },
        "degenerate_by_length_any_first": dict(degenerate_by_length_any_first),
        "case_sensitive": case_sensitive,
        "min_peptide_length": min_peptide_length,
        "ignore_wildcard_peptides": ignore_wildcard_peptides,
        "exclude_single_residue_wildcard_peptides": exclude_single_residue_wildcard_peptides,
        "total_peptides_indexed": indexed_entries
    }


# ============================================================
# Matching
# ============================================================

def format_match(metadata, query_id, query_sequence_original, start_0based):
    length = metadata["length"]

    return {
        "peptide_id": metadata["peptide_id"],
        "peptide_source": metadata["peptide_source"],
        "peptide_sequence": metadata["peptide_sequence"],
        "query_id": query_id,
        "query_sequence": query_sequence_original,
        "match_start": start_0based + 1,
        "match_end": start_0based + length,
        "match_start_0based": start_0based,
        "match_end_exclusive": start_0based + length,
        "matched_query_subsequence": query_sequence_original[start_0based:start_0based + length]
    }


def process_query_record(record):
    global PEPTIDE_INDEX

    query_id, query_sequence_original = record

    if not PEPTIDE_INDEX["case_sensitive"]:
        query_sequence_match = query_sequence_original.upper()
    else:
        query_sequence_match = query_sequence_original

    query_length = len(query_sequence_match)
    matches = []

    exact_by_length = PEPTIDE_INDEX["exact_by_length"]
    degenerate_by_length_and_first = PEPTIDE_INDEX["degenerate_by_length_and_first"]
    degenerate_by_length_any_first = PEPTIDE_INDEX["degenerate_by_length_any_first"]

    all_lengths = set(exact_by_length)
    all_lengths.update(degenerate_by_length_and_first)
    all_lengths.update(degenerate_by_length_any_first)

    for length in sorted(all_lengths):
        if length > query_length:
            continue

        exact_dict = exact_by_length.get(length, {})
        degenerate_first_groups = degenerate_by_length_and_first.get(length, {})
        degenerate_any_first = degenerate_by_length_any_first.get(length, [])

        for start in range(0, query_length - length + 1):
            window = query_sequence_match[start:start + length]

            exact_hits = exact_dict.get(window)
            if exact_hits:
                for metadata in exact_hits:
                    matches.append(
                        format_match(metadata, query_id, query_sequence_original, start)
                    )

            first_residue = window[0]
            candidate_degenerate = degenerate_first_groups.get(first_residue, [])

            if degenerate_any_first:
                candidate_degenerate = candidate_degenerate + degenerate_any_first

            for metadata in candidate_degenerate:
                if metadata["regex"].fullmatch(window):
                    matches.append(
                        format_match(metadata, query_id, query_sequence_original, start)
                    )

    return query_id, matches, None


def chunk_records(records, chunk_size):
    for i in range(0, len(records), chunk_size):
        yield records[i:i + chunk_size]


def process_query_chunk(records):
    chunk_matches = []
    failed = {}

    for record in records:
        query_id = record[0]
        try:
            _, matches, _ = process_query_record(record)
            chunk_matches.extend(matches)
        except Exception as exc:
            failed[query_id] = repr(exc)

    processed_ids = [record[0] for record in records]
    return processed_ids, chunk_matches, failed


# ============================================================
# Longest-overlap post-processing
# ============================================================

def filter_to_longest_overlapping_matches(results_df):
    """
    Remove shorter peptide matches that are fully contained within longer
    peptide matches in the same query sequence.

    Matches at independent locations are retained.
    Multiple annotations with identical coordinates are retained.
    """

    if results_df.empty:
        return results_df

    filtered_groups = []

    for query_id, group in results_df.groupby("query_id", sort=False):
        group = group.copy()

        group["match_length"] = (
            group["match_end_exclusive"] - group["match_start_0based"]
        )

        group = group.sort_values(
            by=["match_length", "match_start_0based", "match_end_exclusive"],
            ascending=[False, True, False]
        )

        kept_rows = []

        for _, row in group.iterrows():
            start = row["match_start_0based"]
            end = row["match_end_exclusive"]
            length = row["match_length"]

            contained_in_longer_match = False

            for kept in kept_rows:
                kept_start = kept["match_start_0based"]
                kept_end = kept["match_end_exclusive"]
                kept_length = kept["match_length"]

                same_coordinates = start == kept_start and end == kept_end

                fully_contained = (
                    start >= kept_start
                    and end <= kept_end
                    and length < kept_length
                )

                if fully_contained and not same_coordinates:
                    contained_in_longer_match = True
                    break

            if not contained_in_longer_match:
                kept_rows.append(row)

        filtered_groups.append(pd.DataFrame(kept_rows))

    filtered_df = pd.concat(filtered_groups, ignore_index=True)

    if "match_length" in filtered_df.columns:
        filtered_df = filtered_df.drop(columns=["match_length"])

    return filtered_df


# ============================================================
# Progress and matching orchestration
# ============================================================

def report_progress(completed, total_queries, start_time):
    percent = (completed / total_queries) * 100 if total_queries else 100
    elapsed = str(datetime.timedelta(seconds=int(time.time() - start_time)))
    print(
        f"[{datetime.datetime.now().strftime('%H:%M:%S')}] "
        f"{completed}/{total_queries} queries completed "
        f"({percent:.2f}%) after {elapsed}"
    )


def find_matching_sequences(
    library_file,
    queries_df,
    case_sensitive=False,
    max_workers=1,
    verbose_every=60,
    progress_path=None,
    min_peptide_length=5,
    chunk_size=10,
    save_every_completed=25,
    ignore_wildcard_peptides=False,
    exclude_single_residue_wildcard_peptides=True
):
    global PEPTIDE_INDEX

    if progress_path:
        os.makedirs(os.path.dirname(progress_path), exist_ok=True)

    peptide_index = build_peptide_index_from_file(
        library_file=library_file,
        case_sensitive=case_sensitive,
        min_peptide_length=min_peptide_length,
        ignore_wildcard_peptides=ignore_wildcard_peptides,
        exclude_single_residue_wildcard_peptides=exclude_single_residue_wildcard_peptides
    )

    PEPTIDE_INDEX = peptide_index

    print(f"Indexed {peptide_index['total_peptides_indexed']} peptides")
    print(f"Minimum peptide length retained: {min_peptide_length}")
    print(f"Ignore wildcard-containing peptides: {ignore_wildcard_peptides}")
    print(
        "Exclude single-residue wildcard peptides: "
        f"{exclude_single_residue_wildcard_peptides}"
    )

    all_results = []
    failed_ids = {}
    processed_ids = set()

    progress_data = load_progress(progress_path)
    if progress_data:
        all_results = progress_data.get("results", [])
        failed_ids = progress_data.get("failed_ids", {})
        processed_ids = set(progress_data.get("processed_ids", []))

    total_queries = len(queries_df)
    completed = len(processed_ids)

    remaining_df = queries_df[~queries_df["Q1"].isin(processed_ids)].copy()
    records = list(remaining_df[["Q1", "Q2"]].itertuples(index=False, name=None))
    chunks = list(chunk_records(records, chunk_size))

    print(f"Resuming from {completed}/{total_queries} completed queries")
    print(f"Remaining query sequences: {len(records)}")
    print(f"Worker processes: {max_workers}")
    print(f"Chunk size: {chunk_size}")

    start_time = time.time()
    last_report_time = time.time()
    last_save_completed = completed

    if max_workers <= 1:
        print("Running in single-process mode to avoid multiprocessing memory duplication.")

        for chunk in chunks:
            processed_chunk_ids, chunk_matches, chunk_failed = process_query_chunk(chunk)

            processed_ids.update(processed_chunk_ids)
            all_results.extend(chunk_matches)
            failed_ids.update(chunk_failed)

            completed = len(processed_ids)
            now = time.time()

            should_save = (
                progress_path
                and completed - last_save_completed >= save_every_completed
            )

            should_report = now - last_report_time >= verbose_every

            if should_save or should_report:
                save_progress_atomic(
                    progress_path,
                    {
                        "results": all_results,
                        "processed_ids": list(processed_ids),
                        "failed_ids": failed_ids
                    }
                )
                last_save_completed = completed

            if should_report:
                report_progress(completed, total_queries, start_time)
                last_report_time = now

    else:
        print(
            "Running in multiprocessing mode. "
            "For very large libraries this may require substantial RAM."
        )

        with ProcessPoolExecutor(
            max_workers=max_workers,
            initializer=init_worker,
            initargs=(peptide_index,)
        ) as executor:

            futures = [executor.submit(process_query_chunk, chunk) for chunk in chunks]

            for future in as_completed(futures):
                processed_chunk_ids, chunk_matches, chunk_failed = future.result()

                processed_ids.update(processed_chunk_ids)
                all_results.extend(chunk_matches)
                failed_ids.update(chunk_failed)

                completed = len(processed_ids)
                now = time.time()

                should_save = (
                    progress_path
                    and completed - last_save_completed >= save_every_completed
                )

                should_report = now - last_report_time >= verbose_every

                if should_save or should_report:
                    save_progress_atomic(
                        progress_path,
                        {
                            "results": all_results,
                            "processed_ids": list(processed_ids),
                            "failed_ids": failed_ids
                        }
                    )
                    last_save_completed = completed

                if should_report:
                    report_progress(completed, total_queries, start_time)
                    last_report_time = now

    if progress_path:
        save_progress_atomic(
            progress_path,
            {
                "results": all_results,
                "processed_ids": list(processed_ids),
                "failed_ids": failed_ids
            }
        )

    print(f"Matching complete: {completed}/{total_queries} queries processed")

    if failed_ids:
        print(f"Warning: {len(failed_ids)} queries failed. See failure log.")

    return pd.DataFrame(all_results), failed_ids


# ============================================================
# Output
# ============================================================

def save_output(results_df, failed_ids, base_filename):
    csv_file = f"{base_filename}.csv"
    txt_file = f"{base_filename}.txt"

    results_df.to_csv(csv_file, index=False)
    results_df.to_csv(txt_file, sep="\t", index=False)

    print("Results saved to:")
    print(f"- {csv_file}")
    print(f"- {txt_file}")

    if failed_ids:
        failed_file = f"{base_filename}_FAILED_QUERIES.csv"
        pd.DataFrame(
            [{"query_id": query_id, "error": error} for query_id, error in failed_ids.items()]
        ).to_csv(failed_file, index=False)

        print("Failed-query log saved to:")
        print(f"- {failed_file}")


def remove_old_progress_files(results_path, fasta_basename, current_progress_path):
    if not os.path.exists(results_path):
        return

    current_name = os.path.basename(current_progress_path)

    for filename in os.listdir(results_path):
        if (
            filename.startswith(fasta_basename)
            and filename.endswith(".progress.pkl")
            and filename != current_name
        ):
            try:
                os.remove(os.path.join(results_path, filename))
            except OSError:
                pass


# ============================================================
# Main
# ============================================================

def main():
    library_file = install_libraries_if_needed()
    fasta_file = select_query_file()

    fasta_basename = os.path.splitext(os.path.basename(fasta_file))[0]
    output_base = os.path.join(results_path, f"{fasta_basename}_RESULTS")

    print("PeptID run configuration")
    print(f"Base directory: {base_path}")
    print(f"Library file: {library_file}")
    print(f"Query FASTA: {fasta_file}")
    print(f"Results directory: {results_path}")
    print(f"Workers: {MAX_WORKERS}")
    print(f"Minimum peptide length: {MIN_PEPTIDE_LENGTH}")
    print(f"Case sensitive: {CASE_SENSITIVE}")
    print(f"Ignore wildcard-containing peptides: {IGNORE_WILDCARD_PEPTIDES}")
    print(
        "Exclude single-residue wildcard peptides: "
        f"{EXCLUDE_SINGLE_RESIDUE_WILDCARD_PEPTIDES}"
    )
    print(
        "Keep only longest overlapping matches: "
        f"{KEEP_ONLY_LONGEST_OVERLAPPING_MATCHES}"
    )

    queries_df = parse_fasta_to_df(fasta_file)
    print(f"Parsed FASTA file: {len(queries_df)} query sequences found")

    if queries_df.empty:
        print("No sequences found in FASTA file. Exiting.")
        return

    fasta_hash = compute_file_hash(fasta_file)
    library_hash = compute_file_hash(library_file)

    combined_hash = hashlib.sha256(
        (
            f"{fasta_hash}_{library_hash}_"
            f"{MIN_PEPTIDE_LENGTH}_{CASE_SENSITIVE}_"
            f"{IGNORE_WILDCARD_PEPTIDES}_"
            f"{EXCLUDE_SINGLE_RESIDUE_WILDCARD_PEPTIDES}_"
            f"{KEEP_ONLY_LONGEST_OVERLAPPING_MATCHES}"
        ).encode()
    ).hexdigest()[:16]

    progress_path = os.path.join(
        results_path,
        f"{fasta_basename}_{combined_hash}.progress.pkl"
    )

    remove_old_progress_files(results_path, fasta_basename, progress_path)

    matching_sequences, failed_ids = find_matching_sequences(
        library_file=library_file,
        queries_df=queries_df,
        case_sensitive=CASE_SENSITIVE,
        max_workers=MAX_WORKERS,
        verbose_every=VERBOSE_EVERY,
        progress_path=progress_path,
        min_peptide_length=MIN_PEPTIDE_LENGTH,
        chunk_size=CHUNK_SIZE,
        save_every_completed=SAVE_EVERY_COMPLETED,
        ignore_wildcard_peptides=IGNORE_WILDCARD_PEPTIDES,
        exclude_single_residue_wildcard_peptides=EXCLUDE_SINGLE_RESIDUE_WILDCARD_PEPTIDES
    )

    if not matching_sequences.empty:
        if KEEP_ONLY_LONGEST_OVERLAPPING_MATCHES:
            before_filter = len(matching_sequences)

            matching_sequences = filter_to_longest_overlapping_matches(
                matching_sequences
            )

            after_filter = len(matching_sequences)

            print(
                f"Removed {before_filter - after_filter} shorter peptide matches "
                f"contained within longer matches from the same query sequence"
            )

        save_output(matching_sequences, failed_ids, output_base)

        if os.path.exists(progress_path):
            os.remove(progress_path)

    else:
        print("No matches to save.")

        if failed_ids:
            save_output(pd.DataFrame(), failed_ids, output_base)

        if os.path.exists(progress_path):
            os.remove(progress_path)


if __name__ == "__main__":
    main()
