Walt Shands
Professor Travis Wheeler 
Transmark benchmark
04 September 2017
Transmark - a benchmark for HMMER


Transmark is a software program for creating artificial ‘background’ DNA sequences, implanting real DNA sequences and artificial decoy DNA sequences (decoy ORFs) in those background sequences, and running homology search programs to measure their sensitivity and selectivity at finding the real DNA sequences.


The benchmark builds background sequences by concatenating nucleotides. A certain number of positive test sequences and a certain number of decoy ORFs are placed into the background sequence of nucleotides. By default Transmark makes a benchmark with background sequence that consists of about 5% decoy ORFs.


Each pseudo-chromosome/background sequence is 100 MB in length by default.  Protein sequences will be ~300 amino acids on average, which is about 1000 bases, so Transmark plants about (0.05 * background sequence length / 1000 ORFs) 5000 decoy sequences in each chromosome.


NOTE many of the benchmark parameters are adjustable, including the length of the background sequence and the percent composition of positive/decoy sequence. This means the total length of the final benchmark chromosome will be the selected length for the background sequence plus the length of the positive sequences plus the length of the decoy sequences, resulting in slightly less percentage decoy sequence then selected by program options (this could be fixed in a future release).


Transmark builds background sequences by generating them using an internally supplied HMM. The program can also use an input database. The benchmark creates a set of training and test DNA sequences from a file of DNA MSAs, and inserts test DNA sequences into the background sequence.


Phmmert and other homology search programs look for test sequences in the artificial chromosomes by using the training sequences to build training multi-sequence alignments (MSA)s. For example in one mode phmmert creates an HMM from a training MSA and uses that to search for positive sequences. In another scenario tblastn creates a consensus sequence from a training MSA and uses that to search the background sequences. 


Transmark starts with the script create-benchmark.sh which accepts an input amino acid MSA file, an input DNA MSA file and a working directory. We used a DNA MSA file constructed from the complete set of  Pfam amino acid MSAs, and the complete amino acid Pfam MSA; which includes sequences from all over the tree of life and which uses non standard translation tables. The DNA MSA file is is then filtered so that only sequences with ORFs that map to the full length of the DNA sequence are kept. This needs to be done because when Transmark uses sequence homology search programs that in turn use the standard translation table and to find ORFs. If non standard translation table sequences are used they can contain stop codons that are located before the end of the ORF and that will confuse the homology search tools. Requiring sequences with ORFs at least as long as the sequence itself ensures the homology search programs can effectively search the training sequences.


Transmark then selects a random set of 7362 alignments from the Pfam DNA MSAs, which is roughly half of the Pfam alignments, to use in creating the training and test sequence sets. This ensures there are enough MSAs to provide a total of 27,531 test sequences for the ten benchmark chromosomes.


The benchmark creates the test and training sequences by reading the input DNA MSA and using single linkage clustering to get training and test sets where every training sequence is less than 60% identical to any test sequence, and where every test sequence is less than 70% identical to any other test sequence. Additionally each training sequence must be at least .01% identical to at least one test sequence.[a]


In order to determine if a training sequence is less than 60% identical to any other training sequence Transmark uses blastn with word size four and evalue 100. Blastn finds alignments for the two training sequences and calculates a percent identity for each alignment. Transmark ensures that only alignments with lengths of at least 90% of the query length are considered when examining percent identity, so that the percent identity of small local alignments are not used. This ensures that we don’t use the percent identity of a small local alignment as representative of the percent identity of the sequences.    


Decoy ORFs are constructed by reading an amino acid MSA file and randomly selecting an MSA and then randomly selecting a sequence. The sequence is then shuffled so it reflects the length and amino acid content of a real ORF but is not a real ORF. Then the amino acid sequence is converted to nucleotides and added to the list of positive test sequences to be embedded in the chromosome. Transmark uses Leucine for an unknown amino acid.


After the list of all positive test sequences and decoy sequences has been created the sequences are inserted into the background nucleotide sequence in randomly selected locations.




[a]not sure if it should be 'at lease .01% identical to every test sequence' see lines 1264 to 1302 in rmark-create.c