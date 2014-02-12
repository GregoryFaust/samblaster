#*samblaster*

Written by Greg Faust (gf4ea@virginia.edu)  
[Ira Hall Lab, University of Virginia](http://faculty.virginia.edu/irahall/)

**Current version:** 0.1.9

Current support for Linux only.

##Summary
*samblaster* is a fast and flexible program for marking duplicates in read-id sorted paired-end SAM files.
It can also optionally output discordant read pairs and/or split read mappings to separate SAM files.

When marking duplicates, *samblaster* will require approximately xG of memory per Y read pairs.

##Installation
*samblaster* is self contained and therefore has no installation dependencies beyond **g++** and **make**.  

**g++ version 4.3 or later is required.**

To manually install *samblaster*, do the following:
~~~~~~~~~~~~~~~~~~
git clone git://github.com/GregoryFaust/samblaster.git
cd samblaster
make
cp samblaster /usr/local/bin/.
~~~~~~~~~~~~~~~~~~

##Usage
See the [SAM File Format Specification](http://samtools.sourceforge.net/SAMv1.pdf) for details about the SAM alignment format.

By default, *samblaster* reads SAM input from **stdin** and writes SAM to **stdout** with duplicates marked. Input SAM file must contain paired end data, contain a sequence header, and be sorted by read ids.
Output SAM file will contain all the alignments in the same order as the input, with duplicates marked with SAM FLAG 0x400.

**COMMON USAGE SCENARIOS:**  

To take input alignments directly from _bwa mem_ and output to _samtools view_ to compress SAM to BAM:
```
bwa mem index samp.r1.fq samp.r2.fq | samblaster | samtools view -Sb - > samp.out.bam
```

To additionally output discordant read pairs and split read alignments:  
```
bwa mem index samp.r1.fq samp.r2.fq | samblaster -e -d samp.disc.sam -s samp.split.sam | samtools view -Sb - > samp.out.bam
```

To pull split reads and discordants read pairs from a pre-existing BAM file with duplicates already marked:  
```
samtools view -h samp.bam | samblaster -a -e -d samp.disc.sam -s samp.split.sam -o /dev/null
```

---
**OPTIONS:**
Default values enclosed in square brackets []
```
Input/Output Options:
-i --input          FILE Input SAM file [stdin].
-o --output         FILE Output SAM file for all input alignments [stdout].
-d --discordantFile FILE Output discordant read pairs to SAM file. By default, no discordant file is produced.
-s --splitterFile   FILE Output split reads to SAM file abiding by paramaters below. By default, no splitter file is produced.

Other Options:
-a --acceptDupMarks      Accept duplicate marks already in input file instead of looking for duplicates in the input.
-e --excludeDups         Exclude reads marked as duplicates from discordant and/or splitter files.
-c --maxSplitCount  INT  Maximum number of split alignments for a read to be output to splitter file.[2]
-m --minNonOverlap  INT  Minimum non-overlaping base pairs between two alignments of a read to be output to splitter file.[20]

-h --help                Print samblaster help to stderr.
   --version             Print samblaster version number to stderr.
```

---
**DUPLICATE IDENTIFICATION:**  
A **duplicate** read pair is defined as a pair that has the same *signature* for each mapped read as a previous read pair in the input SAM file.
The *signature* is comprised of the combination of the sequence name, strand, and starting reference offset of the alignment.
For pairs in which both reads are mapped, both signatures must match.
For pairs in which only one side is mapped (an "orphan"), the signature of the mapped read must match a previously seen orphan.
No doubly unmapped pair will be marked as a duplicate.

---
**DISCORDANT IDENTIFICATION:**  
A **discordant** read pair is one which meets all of the following criteria:
- Both side of the read pair are mapped (neither FLAG 0x4 or 0x8 is set).
- The *properly paired* FLAG is set (0x2).
- Secondary alignments (FLAG 0x100) are never output as discordant, although a discordant pair can have secondary alignments associated with them.
- Duplicate read pairs that meet the above criteria will be output as discordant unless the **-e** option is used.
     
---
**SPLIT READ IDENTIFICATION:**


