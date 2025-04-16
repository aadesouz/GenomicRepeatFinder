# GenomicRepeatFinder
Python-based tandem/interspersed repeat finder for FASTA files via sliding/fixed window analysis. 

Run using `python GenomicRepeatFinder.py`

The python script used to identify repeating patterns in genomic sequences requires the following modules:  Biopython’s SeqIO, `re` , `sys` , and `string`. A FASTA file with the genome of interest is required. 

**User defined parameters:**

-       BIN_SIZE = Groups hits based on their location in the genome for graphing.

-       FORWARD_SEQ = The sequence of interest (e.g., `TTTAGGG`).

-       MIN_REPEATS = Minimum number of sequence repeats to qualify. 

-       FASTA_FILE = Path to the FASTA file.

-       WINDOW_SIZE = Length of the analysis window 

-       STEP_SIZE = Step size for adjusting sliding/overlapping window analysis.

###### Functions:
-   `prompt_user()`
	Gathers and verifies user defined parameters

-   `process_all_records(fasta_records, FORWARD_SEQ, MIN_REPEATS, STEP_SIZE, WINDOW_SIZE, define_pattern, circle_match)`
	Loops through each sequence entry of the FASTA file, calling helper functions to complete data collection.

-   `complement(sequence)`
	Returns the complementary DNA sequence (reverse complement) for a given sequence. It uses a dictionary to map bases (A↔T, C↔G).

-   `define_pattern(pattern)`
	Generates a list of regular expression patterns for detecting repeats of a forward sequence. The pattern is designed to allow for variable starting positions and a specified minimum repeat count.

-   `circle_match(pattern, fasta, pos, MIN_REPEATS)`
	Finds matches for a given pattern in the sequence using regular expressions. Returns the number of matching sequences, provided they meet the minimum repeat. 

-   `get_bin(position)`
	Determines which bin a given position in the sequence belongs to by dividing the position by BIN_SIZE.

###### Input collection & FASTA Parsing: 
The main function uses `prompt_user()` to collect parameters. `The SeqIO.parse()` from Biopython is used to parse the fasta file into fasta_records. For each sequence record in fasta_records, we search and document matches using `process_all_records()`. 
###### Pattern Generation & Reverse Complement: 
From the defined sequence of interest, FORWARD_SEQ, the reverse complement sequence is generated using `complement(sequence)`. Using `define_pattern()`, the permutations and their complement are created according to the desired minimum repeats. 
For example, generate and search for each of these 7 versions and their reverse complements. Each version must contain 3 or more repeats of the 7-bp sequence (because 21 bp is required for small RNA to bind). Search all contiguous repeats (including repeats of 3, 4, 5, 6, ..., etc). 

```
              5'-TTTAGGG-3' 5'-TTTAGGG-3' 5'-TTTAGGG-3'
         G-3' 5'-TTTAGGG-3' 5'-TTTAGGG-3' 5'-TTTAGG
        GG-3' 5'-TTTAGGG-3' 5'-TTTAGGG-3' 5'-TTTAG
       GGG-3' 5'-TTTAGGG-3' 5'-TTTAGGG-3' 5'-TTTA
      AGGG-3' 5'-TTTAGGG-3' 5'-TTTAGGG-3' 5'-TTT
     TAGGG-3' 5'-TTTAGGG-3' 5'-TTTAGGG-3' 5'-TT
    TTAGGG-3' 5'-TTTAGGG-3' 5'-TTTAGGG-3' 5'-T

5'-CCCTAAA-3' 5'-CCCTAAA-3' 5'-CCCTAAA-3'
    CCTAAA-3' 5'-CCCTAAA-3' 5'-CCCTAAA-3' 5'-C
     CTAAA-3' 5'-CCCTAAA-3' 5'-CCCTAAA-3' 5'-CC
      TAAA-3' 5'-CCCTAAA-3' 5'-CCCTAAA-3' 5'-CCC
       AAA-3' 5'-CCCTAAA-3' 5'-CCCTAAA-3' 5'-CCCT
        AA-3' 5'-CCCTAAA-3' 5'-CCCTAAA-3' 5'-CCCTA
         A-3' 5'-CCCTAAA-3' 5'-CCCTAAA-3' 5'-CCCTAA
```

###### Pattern Matching & Binning Matches: 
A dictionary, observations, is initiated to track forward and reverse hits in each bin. The script slides across the sequence with STEP_SIZE, extracting from windows of length, WINDOW_SIZE, making sliding window analysis possible. In each window, the script matches forward and reverse patterns using the `circle_match()` function, which uses `re.finditer()`. Matches can be found starting from any position within the window.
A dictionary, observations, is initiated to track forward and reverse hits in each bin. The script slides across the sequence with STEP_SIZE, extracting from windows of length, WINDOW_SIZE, making sliding window analysis possible. In each window, the script matches forward and reverse patterns using the `circle_match()` function, which uses `re.finditer()`. Matches can be found starting from any position within the window.



Many thanks to Dr. Li, Dawei in the completion of this program. 
