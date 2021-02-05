# Smith-Waterman
Implementation of the Smith–Waterman algorithm

### Basic Usage

```sh
  $ python smith_waterman.py --help
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
