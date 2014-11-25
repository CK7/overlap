##overlap

overlap is a program that identifies overlaps between long reads that can be used for the assembly of the reads (using the software Lola, see 
separate repository) or for other purposes. Program input consists of a blast file with the input sequences aligned. The program will use hits in the blast file as 
seeds and will try to elongate them using different versions of the Smith-Waterman algorithm specifically designed to identify end-end and containment 
overlaps. 

### Running overlap

In order to run overlap you need to

1. Run blast on the sequences. You should use the -m 8 option. We recommend that you also use -F F -r 1 -q -5 -e 1e-30. The -e parameter will prevent blast from inflating the output file. 
2. Run overlap

For running overlap use the following command:

```
$ <path-to-Lola-directory>/overlap
```

You can get a description of the parameters taken by overlap by running this command without parameters:

```
$ ./overlap

./overlap v1.01

Usage: overlap [-os <min-overlap-size>] [-G <gap-penalty>] [-r <match-reward>] [-d <seed-dize>]
               [-q <mismatch-penalty>] [-p <% identity>] <seq-file> <m8-blast-file>

Where
 -os      | minimum overlap required for determining connection (default: 500 bp)
 -d       | minimum alignment size required for a couple to be considered in the m8 blast file (default: 300)
 -p       | % identity threshold for a connection (default: 99%)
 -G       | Cost to open a gap (defualt: -5)
 -r       | Reward for a nucleotide match (default: 1)
 -q       | Penalty for a nucleotide mismatch (default: -3)
 m8-blast-file is a self-blast report (recommended aruments: -F F -m 8 -r 1 -q -3)
 seq-file is a FASTA file of the analyzed sequences
```

### overlap output
Output is written to a file named **\<seq-file-without-postfix\>.overlap.txt** and contains a report for all overlaps detected. File consists of line with 
two types of format:

```
<TYPE1>	<seq1>	<side1>	<seq2>	<side2>	<start1>	<end1>	<start2>	<end2>	<size1>	<size2>	<% identity>
```
where
* TYPE1 can be either CONNECTED (the two sequence have an end to end overlap) or SHARED (at least one of the sequences has a non-end region overlapping with 
the other sequence or the two ends overlap in a way that does not enable their assembly)
* side1 and side2 can be either 3 or 5 (for CONNECTED) or also Middle (for SHARED)
* start1 and end1 are start and end coordinates of the overlapping region for seq1
* start2 and end2 are coordinates for overlapping region of seq2. start2 can be smaller than end2 if reverse complement of seq2 is aligned
* size1 and size2 are the lengths of seq1 and seq2, repectively
* % identity is % identity for the alignment

The second format is
```
<TYPE2>	<seq1>	<seq2>	<start1>	<end1>	<start2>	<end2>	<size1>	<size2>	<% identity>
```

where
* TYPE2 can be either IDENTICAL (the two sequences align throughout their whole length except maybe for a few bps) or CONTAINS (seq1 contains seq2)
* Other fields are the same as above

### example

Refer to the directory example in the Lola repository for an example of running overlap. File example/RBG-1.reads.overlap.txt is the output of running 
overlap on the file RBG-1.reads.fna. See example/run for command line used to generate this file. 
