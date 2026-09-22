PeptID

PeptID is a Python tool for identifying known peptide sequences within protein FASTA files. It compares query sequences against the PeptID combined peptide library, records the coordinates and source annotation of every match, and exports the results in both CSV and tab-delimited formats.

The program is designed to be portable: it creates its own working directories and downloads the required reference libraries automatically from the latest PeptID GitHub release.

Key features
Exact peptide matching against protein or proteome FASTA files

Optional support for X wildcard residues in library peptides

Configurable minimum peptide length and case sensitivity

Removal of shorter matches contained entirely within longer matches

Automatic downloading and extraction of the PeptID reference libraries

Automatic detection of a query file in the Queries directory

Progress checkpointing for long analyses

Optional multiprocessing

CSV and tab-delimited result files

Windows, Linux, and macOS-compatible paths

Requirements
Python 3.9 or newer is recommended

pandas

An internet connection when the libraries are installed for the first time

Sufficient disk space for both the compressed and extracted reference libraries

Install the required Python package with:

python -m pip install pandas
On some systems, the Python command may be python3 instead of python.

Downloading PeptID
Download as a ZIP
Open the PeptID GitHub repository.

Select Code and then Download ZIP.

Extract the downloaded ZIP to a location where you have write permission.

Clone with Git
git clone https://github.com/SRF-Boltz/PeptID.git
cd PeptID
Quick start
Place PeptID.py in its own directory.

Create a directory named Queries beside the script.

Put one protein FASTA file in Queries.

Open a terminal or command prompt in the PeptID directory.

Run:

python PeptID.py
For example:

PeptID/
├── PeptID.py
└── Queries/
    └── example_proteome.fasta
PeptID will create the remaining directories, download and extract the reference libraries if necessary, analyse the query sequences, and write the results to Results.

You can also run PeptID once before adding a query. It will create the required directories and install the libraries, then stop with a message asking you to place a query file in Queries.

Automatic library installation
On first use, PeptID downloads these assets from the latest published release of SRF-Boltz/PeptID:

Function_Library.zip

Source_Library.zip

Combined_Library.zip

All three archives are saved and extracted in Libraries. PeptID uses the extracted Combined_Library.txt as its matching database.

After installation, a typical PeptID directory looks like:

PeptID/
├── PeptID.py
├── Libraries/
│   ├── Function_Library.zip
│   ├── Source_Library.zip
│   ├── Combined_Library.zip
│   ├── Combined_Library.txt
│   └── [other extracted library content]
├── Queries/
│   └── example_proteome.fasta
└── Results/
The precise extracted layout may vary if an archive contains a subdirectory. PeptID searches the Libraries directory recursively for Combined_Library.txt.

On later runs, PeptID reuses the installed libraries and skips the download. To perform a clean reinstallation, remove the Libraries directory and run the program again.

Query input
PeptID accepts one FASTA-like file in Queries with any of these extensions:

.fasta  .fa  .faa  .fas  .fna  .txt
The file must contain FASTA headers beginning with > followed by one or more sequence lines:

>protein_1
MSTNPKPQRKTKRNTNRRPQDVKFPGGGQIVGGVYLLPRRG
>protein_2
MALWMRLLPLLALLALWGPGPGAGSLQPLALEGSLQKRGIVE
By default, there must be exactly one supported query file in Queries. If PeptID finds more than one, it stops rather than guessing which file to analyse.

To select a particular file while keeping several files in the directory, edit this setting near the top of PeptID.py:

QUERY_FILENAME = "example_proteome.fasta"
Leave it as None for automatic detection:

QUERY_FILENAME = None
Results
For a query named example_proteome.fasta, PeptID writes:

Results/
├── example_proteome_RESULTS.csv
└── example_proteome_RESULTS.txt
The .txt file is tab-delimited. Both output files contain the same results.

Column	Description
peptide_id	Identifier from the PeptID library
peptide_source	Source annotation from the library
peptide_sequence	Peptide sequence stored in the library
query_id	FASTA header of the matched query sequence
query_sequence	Complete query sequence
match_start	One-based inclusive start position
match_end	One-based inclusive end position
match_start_0based	Zero-based start position
match_end_exclusive	Zero-based exclusive end position
matched_query_subsequence	Sequence recovered from the matched region
If one or more query sequences fail during processing, PeptID also creates:

example_proteome_RESULTS_FAILED_QUERIES.csv
Matching behaviour
The default configuration is:

MAX_WORKERS = 1
CHUNK_SIZE = 10
MIN_PEPTIDE_LENGTH = 5
CASE_SENSITIVE = False
IGNORE_WILDCARD_PEPTIDES = True
EXCLUDE_SINGLE_RESIDUE_WILDCARD_PEPTIDES = True
KEEP_ONLY_LONGEST_OVERLAPPING_MATCHES = True
Peptide length
Library entries shorter than MIN_PEPTIDE_LENGTH are excluded before matching.

Case sensitivity
Matching is case-insensitive by default. Set CASE_SENSITIVE = True if letter case should be treated as meaningful.

Wildcard residues
PeptID can interpret X in a library peptide as a single-position wildcard. Wildcard-containing peptides are ignored under the default setting:

IGNORE_WILDCARD_PEPTIDES = True
To enable them, change this to False. When wildcard matching is enabled, the separate EXCLUDE_SINGLE_RESIDUE_WILDCARD_PEPTIDES option prevents highly nonspecific peptides containing only one defined residue from being used.

Overlapping matches
With KEEP_ONLY_LONGEST_OVERLAPPING_MATCHES = True, a shorter match is removed when it is entirely contained within a longer match in the same query sequence. Matches at independent positions are retained, as are multiple library annotations with identical coordinates.

Interrupted analyses and progress files
PeptID periodically saves a checkpoint in Results using a .progress.pkl file. If an analysis is interrupted, rerunning PeptID with the same query, library, and matching settings resumes from the saved checkpoint.

The checkpoint filename incorporates hashes of the query and combined library together with the relevant matching settings. This prevents an incompatible checkpoint from being reused accidentally. The checkpoint is removed after a successful run.

Multiprocessing and memory use
PeptID defaults to:

MAX_WORKERS = 1
This avoids duplicating a large peptide index across worker processes. Systems with sufficient RAM can use more processes by increasing MAX_WORKERS, but memory use may rise substantially.

Troubleshooting
No module named 'pandas'
Install pandas using:

python -m pip install pandas
No query FASTA was found
Place one supported FASTA-like file in the Queries directory and rerun PeptID.

Multiple query files were found
Remove the additional query files or set QUERY_FILENAME to the exact filename you want to analyse.

A library download returns HTTP 404
Confirm that:

the PeptID repository is public;

at least one GitHub release has been published; and

the release contains all three assets with the exact filenames listed above.

Draft releases and ordinary repository files are not used by the automatic downloader.

A download was interrupted
Run PeptID again. Incomplete .part files are discarded automatically, while complete valid ZIP archives are reused.

Permission errors
Move the PeptID directory to a location where your user account can create and modify files. Avoid read-only or administrator-controlled installation directories.

Reproducibility
PeptID downloads library assets from the latest published release when a new installation is performed. Record the PeptID release tag and library version used for each analysis. Retain the generated result files and relevant configuration values with the associated research record.

Reporting problems
Please report reproducible problems through the repository's GitHub Issues page. Include:

your operating system;
your Python version;
the complete error message;
the relevant PeptID configuration values; and
a minimal example query where possible.
Do not upload confidential or unpublished biological sequence data to a public issue.
