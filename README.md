# Primer-Amplification

## Description
The program uses a file `sequence.fasta`, containing the sequence of nucleic acids, and `Introns.csv`, having positions of intronic sequence, to generate the following : 
-  All possible reading frames for the sequence given.
-  Reading frame responsible for getting the protein sequence
-  Protein sequence generated from the above reading frame
-  Using the rules given in the Rules_primer.doc generate 3 sets of primers for amplification of sequence to be cloned into Pichia Pastoris.
-  Determine what nucleotide sequence will you use to clone given sequence in the following hosts for large therapeutic applications.
		a) E.Coli.
		b) Pichia Pastoris.
		c) HEK293 Cell line.

## Running the program
- Run the command `python3 exec.py` to generate 5 files each of which corresponds to the description points as mentioned.
