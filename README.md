#samblaster

**Current version:** 0.1.9

Current support for Linux only

Created by Greg Faust (gf4ea@virginia.edu)  
[Ira Hall Lab, University of Virginia](http://faculty.virginia.edu/irahall/)

##Summary
samblaster is a program for quicky processing read-id sorted paired-end SAM files.  It can mark duplicate read pairs, and optionally simultaneously output discordant read pairs and/or split read mappings to separate sam files.

##Installation
samblaster is self contained and therefore has no installation dependencies beyond **g++** and **make**.  

**g++ version 4.3 or later is required.**

To manually install, do the following:
~~~~~~~~~~~~~~~~~~
git clone git://github.com/GregoryFaust/samblaster.git
cd samblaster
make
cp samblaster /usr/local/bin/.
~~~~~~~~~~~~~~~~~~

##Usage
By default, samblaster is a tool that takes SAM input from stdin and write SAM to stdout with duplicates marked (flag 0x400 set).

samblaster -h will print the following help message:

Tool to mark duplicates and optionally output split reads and/or discordant pairs.
Input sam file must contain paired end data, contain sequence header and be sorted by read ids.
Output will be all alignments in the same order as input, with duplicates marked with FLAG 0x400.

Usage:
For use with as a post process on an aligner (eg. bwa mem):
     bwa mem index samp.r1.fq samp.r2.fq | samblaster [-e] [-d samp.disc.sam] [-s samp.split.sam] | samtools view -Sb - > samp.out.bam
For use with a pre-existing bam file to pull split reads and/or discordants:
     samtools view -h samp.bam | samblaster [-e] [-d samp.disc.sam] [-s samp.split.sam] -o /dev/null

Input/Output Options:
-i --input          FILE Input sam file [stdin].
-o --output         FILE Output sam file for all input alignments [stdout].
-d --discordantFile FILE Output discordant read pairs to this file. By default, no discordant file is produced.
-s --splitterFile   FILE Output split reads to this file abiding by paramters below. By default, no splitter file is produced.

Other Options:
-a --acceptDupMarks      Accept duplicate marks already in input file instead of looking for duplicates in the input.
-e --excludeDups         Exclude reads marked as duplicates from discordant and/or splitter file.
-c --maxSplitCount  INT  Maximum number of split alignments for a read to be included in splitter file.[2]
-m --minNonOverlap  INT  Minimum non-overlaping base pairs between two alignments for a read to be included in splitter file. [20]

### Duplicate Identification:
A **duplicate** read pair is defined as a pair that has the same seq/strand/starting-reference-offset *signature* for each mapped read as a previously seen read pair.
For pairs in which both reads are mapped, both signatures much match.
For pairs in which only one side is mapped (an "orphan"), the signature of the mapped read must match a previously seen orphan.
No doubly unmapped pair will be marked as a duplicate.

### Discordant Identification:
A **discordant** read pair is one which meets all of the following criteria:
- Both side of the read pair are mapped (neither flag 0x4 or 0x8 is set).
- The **properly paired** flag is set (0x02).
- Secondary alignments (flag 0x100) are never output as a discordant, although a discordant pair can have secondary alignments associated with them.

     
### Split Read Identification:


