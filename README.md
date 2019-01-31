# Needleman-Wunsch-Algortihm

The Needleman–Wunsch algorithm is an algorithm used in bioinformatics to align protein or nucleotide sequences.

## Usage

Needleman.py contains the main method. The program takes command-line arguments with
the following options.

```
python Needleman.py sequences_file_path substitution_matrix_file_path gap_penalty
```

## Input

### Sequences
The sequences file must be formatted in that way:

```
>sequence 1 name
SEQUENCE1
>sequence 2 name
SEQUENCE2
```

### Substitution matrix
The substitution matrix describes the rate at which one character in a sequence changes to other character states over time.
Other substitution matrices exist, here the Pam250 substitution matrix is provided.

### Gap penalty

a Gap penalty is a method of scoring alignments of two or more sequences. When aligning sequences, introducing a gaps in the sequences can allow an alignment algorithm to match more terms than a gap-less alignment can.

## Example

```
python Needleman.py seq.txt pam250.tab -6
```

