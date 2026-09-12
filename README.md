# Smith-Waterman

Implementation of the **Smith-Waterman algorithm** for local alignment of two biological sequences.

The program accepts **DNA/nucleotide** or **amino acid** sequences in a FASTA file. For protein sequences, a BLOSUM or PAM substitution matrix can be selected. For nucleotide sequences, the `match` and `mismatch` values can be defined.

## Requirements

- Python 3
- Dash
- Dash Bio

Install the required Python packages with:

```bash
pip install dash
pip install dash-bio
```

## Usage

```bash
python smith_waterman.py -t TYPE -f FILE [OPTIONS]
```

```sh
python smith_waterman.py --help
  usage: smith_waterman.py [-h] -t {nt,aa} [-sm {BLOSUM45,BLOSUM50,BLOSUM62,BLOSUM80,BLOSUM90,PAM30,PAM70,PAM250}] -f FILE [-m MATCH] [-mi MISMATCH_PENALTY] [-gap GAP_PENALTY] [-o FOLDER] [--version]

  Implementation of the Smith–Waterman algorithm

  optional arguments:
    -h, --help            show this help message and exit
    -t {nt,aa}, --type {nt,aa}
                          nt: Nucleotide sequence | aa: Amino acid sequence
    -sm {BLOSUM45,BLOSUM50,BLOSUM62,BLOSUM80,BLOSUM90,PAM30,PAM70,PAM250}, --substitution_matrix {BLOSUM45,BLOSUM50,BLOSUM62,BLOSUM80,BLOSUM90,PAM30,PAM70,PAM250}
                          Substitution Matrix type (Only for amino acid sequence) [default: BLOSUM62].
    -f FILE, --fasta FILE
                          Fasta file
    -m MATCH, --match MATCH
                          Match value (Only for nucleotide sequence) [default: 1].
    -mi MISMATCH_PENALTY, --mismatch_penalty MISMATCH_PENALTY
                          Mismatch penalty value (Only for nucleotide sequence) [default: 0].
    -gap GAP_PENALTY, --gap_penalty GAP_PENALTY
                          Gap penalty value [default: 0].
    -o FOLDER, --output FOLDER
                          Output folder
    --version             show program's version number and exit

  Examples of alignment:
    For amino acid sequences
      python smith_waterman.py -t aa -f sequences.fa -gap -1

    For nucleotide sequences
      python smith_waterman.py -t nt -f sequences.fa -m 2 -mi -1 -gap -2

  Thank you!
```

### Parameters

| Parameter | Description | Values |
|---|---|---|
| `-t`, `--type` | Sequence type | `nt` or `aa` |
| `-f`, `--fasta` | FASTA file containing the sequences | FASTA file |
| `-sm`, `--substitution_matrix` | Substitution matrix for amino acid sequences | BLOSUM45, BLOSUM50, BLOSUM62, BLOSUM80, BLOSUM90, PAM30, PAM70, or PAM250 |
| `-m`, `--match` | Score for a nucleotide match | integer |
| `-mi`, `--mismatch_penalty` | Penalty for a nucleotide mismatch | integer |
| `-gap`, `--gap_penalty` | Gap penalty | integer |
| `-o`, `--output` | Output directory | directory |
| `--version` | Show the program version | — |

The default substitution matrix for amino acid sequences is **BLOSUM62**. For nucleotide sequences, the default `match` value is `1` and the default `mismatch` value is `0`.

## Examples

### DNA sequences

```bash
python smith_waterman.py -t nt -f file.fa -m 2 -mi -1 -gap -3 -o out_align
```

### Amino acid sequences

```bash
python smith_waterman.py -t aa -f file.fa -sm BLOSUM62 -gap -2 -o out_align
```

## Input FASTA

The FASTA file should contain two sequences. For example:

### DNA sequences
```text
>sequence1
AATTTACGCGGCATTATAGATACAATCGTGTCT
>sequence2
GCAATTGGCCGGAATTTAATTGATACAGCGC
```

### Amino acid sequences
```text
>protein1
MKTIIALSYIFCLVFADYKDDDDK
>protein2
MKTIIALSYIFCLVFADYKDEDDK
```

The program reads the first two sequences from the FASTA file.

## Output

The program reports:

- parameters used;
- alignment score;
- sequence alignment.

The output directory also contains:

- `alignment_matrix.txt` — scoring matrix generated during the alignment;
- `log_smith_waterman_YYYYMMDD.log` — execution log.

The resulting alignment is **local**, so it can correspond to only part of the input sequences.
