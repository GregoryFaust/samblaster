#*samblaster*

**Written by:** Greg Faust (gf4ea@virginia.edu)  
[Ira Hall Lab, University of Virginia](http://faculty.virginia.edu/irahall/)

**Please cite:**  
[Faust, G.G. and Hall, I.M., “*SAMBLASTER*: fast duplicate marking and structural variant read extraction,” *Bioinformatics* Sept. 2014; **30**(17): 2503-2505.](http://bioinformatics.oxfordjournals.org/content/30/17/2503)

**Also see:** [SAMBLASTER_Supplemental.pdf] (https://www.dropbox.com/s/s7vvdtf2gmytvay/SAMBLASTER_Supplemental.pdf?dl=0) for additonal discussion and statistics about the duplicates marked by *samblaster* vs. *Picard* using the NA12878 sample dataset.  
Click the preceeding link or download the file from this repository.

---

**Current version:** 0.1.24

Support for Linux and OSX (Version 10.7 or higher).

##Summary
*samblaster* is a fast and flexible program for marking duplicates in __read-id grouped<sup>1</sup>__ paired-end SAM files.
It can also optionally output discordant read pairs and/or split read mappings to separate SAM files, and/or unmapped/clipped reads to a separate FASTQ file.
When marking duplicates, *samblaster* will require approximately 20MB of memory per 1M read pairs.

##Installation
*samblaster* is self contained and therefore has no installation dependencies beyond **g++** and **make**.  

*samblaster* can be downloaded from the **_releases_** tab or manually downloaded via *git clone*.  Afterwards, simply use *make* and copy *samblaster* to a directory in your path.  For example:
~~~~~~~~~~~~~~~~~~
git clone git://github.com/GregoryFaust/samblaster.git
cd samblaster
make
cp samblaster /usr/local/bin/.
~~~~~~~~~~~~~~~~~~

##Usage
See the [SAM File Format Specification](http://samtools.sourceforge.net/SAMv1.pdf) for details about the SAM alignment format.

By default, *samblaster* reads SAM input from **stdin** and writes SAM to **stdout**. Input SAM files usually contain paired end data (see [Duplicate Identification](#DupIdentification) below), must contain a sequence header, and must be __read-id grouped<sup>1<sup>__.
By default, the output SAM file will contain all the alignments in the same order as the input, with duplicates marked with SAM FLAG 0x400.  The **--removeDups** option will instead remove duplicate alignments from the output file.

__<sup>1</sup>A read-id grouped__ SAM file is one in which all alignments for a read-id (QNAME) are grouped together in adjacent lines.
Aligners naturally produce such files.
They can also be created by sorting a SAM file by read-id. 
But as shown below, sorting the input to *samblaster* by read-id is not required if the alignments are already grouped.

**COMMON USAGE SCENARIOS:**  

To take input alignments directly from _bwa mem_ and output to _samtools view_ to compress SAM to BAM:
```
bwa mem <idxbase> samp.r1.fq samp.r2.fq | samblaster | samtools view -Sb - > samp.out.bam
```

When using the *bwa mem* **-M** option, also use the *samblaster* **-M** option:
```
bwa mem -M <idxbase> samp.r1.fq samp.r2.fq | samblaster -M | samtools view -Sb - > samp.out.bam
```

To additionally output discordant read pairs and split read alignments:  
```
bwa mem <idxbase> samp.r1.fq samp.r2.fq | samblaster -e -d samp.disc.sam -s samp.split.sam | samtools view -Sb - > samp.out.bam
```

To pull split reads and discordants read pairs from a pre-existing BAM file with duplicates already marked:  
```
samtools view -h samp.bam | samblaster -a -e -d samp.disc.sam -s samp.split.sam -o /dev/null
```

---
**OPTIONS:**
Default values enclosed in square brackets []
<pre>
Input/Output Options:
-i --input           FILE Input sam file [stdin].
-o --output          FILE Output sam file for all input alignments [stdout].
-d --discordantFile  FILE Output discordant read pairs to this file. [no discordant file output]
-s --splitterFile    FILE Output split reads to this file abiding by paramaters below. [no splitter file output]
-u --unmappedFile    FILE Output unmapped/clipped reads as FASTQ to this file abiding by parameters below. [no unmapped file output].
                          Requires soft clipping in input file.  Will output FASTQ if QUAL information available, otherwise FASTA.

Other Options:
-a --acceptDupMarks       Accept duplicate marks already in input file instead of looking for duplicates in the input.
-e --excludeDups          Exclude reads marked as duplicates from discordant, splitter, and/or unmapped file.
-r --removeDups           Remove duplicates reads from all output files. (Implies --excludeDups).
   --addMateTags          Add MC and MQ tags to all output paired-end SAM lines.
   --ignoreUnmated        Suppress abort on unmated alignments. Use only when sure input is read-id grouped and alignments have been filtered.
                          <b>--ignoreUnmated is not recommended for general use. It disables checks that detect incorrectly sorted input.</b>
-M                        Compatibility mode (details below); both FLAG 0x100 and 0x800 denote supplemental (chimeric). Similar to <i>bwa mem</i> <b>-M</b> option.
   --maxSplitCount    INT Maximum number of split alignments for a read to be included in splitter file. [2]
   --maxUnmappedBases INT Maximum number of un-aligned bases between two alignments to be included in splitter file. [50]
   --minIndelSize     INT Minimum structural variant feature size for split alignments to be included in splitter file. [50]
   --minNonOverlap    INT Minimum non-overlaping base pairs between two alignments for a read to be included in splitter file. [20]
   --minClipSize      INT Minumum number of bases a mapped read must be clipped to be included in unmapped file. [20]

-h --help                 Print samblaster help to stderr.
-q --quiet                Output fewer statistics.
   --version              Print samblaster version number to stderr.
</pre>

---
**ALIGNMENT TYPE DEFINITIONS:<a name="Definitions"></a>**  
Below, we will use the following definitions for alignment types.
Starting with *samblaster* release 0.1.22, these definitions are affected by the use of the **-M** option.
By default, *samblaster* will use the current definitions of alignment types as specified in the [SAM Specification](http://samtools.sourceforge.net/SAMv1.pdf).
Namely, alignments marked with FLAG 0x100 are considered *secondary*, while those marked with FLAG 0x800 are considered *supplemental*.
If the **-M** option is specified, alignments marked with either FLAG 0x100 or 0x800 are considered *supplemental*, and no alignments are considered *secondary*.
A *primary* alignment is always one that is neither *secondary* nor *supplemental*.
Only *primary* and *supplemental* alignments are used to find chimeric (split-read) mappings.
The **-M** flag is used for backward compatibility with older SAM/BAM files in which "chimeric" alignments were marked with FLAG 0x100, and should also be used with output from more recent runs of *bwa mem* using its **-M** option.

---
**DUPLICATE IDENTIFICATION:<a name="DupIdentification"></a>**  
A **duplicate** read pair is defined as a pair that has the same *signature* for each mapped read as a previous read pair in the input SAM file.  The *signature* is comprised of the combination of the sequence name, strand, and the reference offset where the 5' end of the read would fall if the read were fully aligned (not clipped) at its 5' end.  The 5' aligned reference position is calculated using a combination of the POS field, the strand, and the CIGAR string.  This definition of *signature* matches that used by *Picard MarkDuplicates*.

1. For pairs in which both reads are mapped, both signatures must match.
2. For pairs in which only one side is mapped (an "orphan"), the signature of the mapped read must match a previously seen orphan. In an orphan pair, the unmapped read need not appear in the input file. In addition, mapped non-paired single read alignments will be treated the same as an orphan pair with a missing unmapped read.
3. No doubly unmapped pair will be marked as a duplicate.
4. Any *secondary* or *supplemental* alignment associated with a duplicate *primary* alignment will also be marked as a duplicate.

---
**DISCORDANT READ PAIR IDENTIFICATION:**  
A **discordant** read pair is one which meets all of the following criteria:

1. Both side of the read pair are mapped (neither FLAG 0x4 or 0x8 is set).
2. The *properly paired* FLAG (0x2) is not set.
3. *Secondary* or *supplemental* alignments are never output as discordant, although a discordant read pair can have such alignments associated with them.
4. Duplicate read pairs that meet the above criteria will be output as discordant unless the **-e** option is used.
     
---
**SPLIT READ IDENTIFICATION:**  
**Split Read** alignments are derived from a single read when one portion of the read aligns to a different region of the reference genome than another portion of the read.  Such pairs of alignments often define a structural variant (SV) breakpoint, and are therefore useful input to SV detection algorithms such as [LUMPY](https://github.com/arq5x/lumpy-sv/).  *samblaster* uses the following strategy to identify split reads alignments.

1. Identify reads that have between two and **--maxSplitCount** *primary* and *supplemental* alignments.
2. Sort these alignments by their strand-normalized position along the read.
3. Two alignments are output as splitters if they are adjacent on the read, and meet these criteria:
    - each covers at least **--minNonOverlap** base pairs of the read that the other does not.
    - the two alignments map to different reference sequences and/or strands.
    - the two alignments map to the same sequence and strand, and represent a SV that is at least **--minIndelSize** in length, and have at most **--maxUnmappedBases** of un-aligned base pairs between them.
4. Split read alignments that are part of a duplicate read will be output unless the **-e** option is used.

---
**UNMAPPED/CLIPPED READ IDENTIFICATION:** 
An **unmapped** or **clipped** read is a *primary* alignment that is unaligned over all or part of its length respectively.  The lack of a full alignment may be caused by a SV breakpoint that falls within the read.  Therefore, *samblaster* will optionally output such reads to a FASTQ file for re-alignment by a tool, such as [YAHA](https://github.com/GregoryFaust/yaha/), geared toward finding split-read mappings.  *samblaster* applies the following strategy to identify and output unmapped/clipped reads:

1. An **unmapped** read has the *unmapped read* FLAG set (0x4).
2. A **clipped** read is a mapped read with a CIGAR string that begins or ends with at least **--minClipSize** unaligned bases (CIGAR code S and/or H), and is not from a read that has one or more *supplemental* alignments.
3. In order for *samblaster* to output the entire sequence for clipped reads, the input SAM file must have soft clipped primary alignments.
4. *samblaster* will output unmapped/clipped reads into a FASTQ file if QUAL information is available in the input file, and a FASTA file if not.
5. Unmapped/clipped reads that are part of a duplicate read pair will be output unless the **-e** option is used.
