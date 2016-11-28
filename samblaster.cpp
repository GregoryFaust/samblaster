/* -*- mode: C++ ; indent-tabs-mode: nil ; c-file-style: "stroustrup" -*-

    Project: samblaster
             Fast mark duplicates in read-ID grouped SAM file.
             Also, optionally pull discordants, splitters, and/or unmappend/clipped reads.
    Author:  Greg Faust (gf4ea@virginia.edu)
    Date:    October 2013

    File:    samblaster.cpp  code file for the main routine and most of the other code.

    License Information:

    Copyright 2013-2016 Gregory G. Faust

    Licensed under the MIT license (the "License");
    You may not use this file except in compliance with the License.
    You may obtain a copy of the License at http://opensource.org/licenses/MIT

*/

// This define is needed for portable definition of PRIu64
#define __STDC_FORMAT_MACROS

#include <stdlib.h>
#include <inttypes.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <map>
#include "sbhash.h"

// Rename common integer types.
// I like having these shorter name.
typedef uint64_t UINT64;
typedef uint32_t UINT32;
typedef int32_t   INT32;

// Some helper routines.

// mempcpy is a GNU extension and not available everywhere.
#ifndef _GNU_SOURCE
inline void *mempcpy(void *dest, const void *src, size_t n)
{
    return (char*) memcpy(dest, src, n) + n;
}
#endif

inline bool streq(char * s1, const char * s2) __attribute__((always_inline));
inline bool streq(char * s1, const char * s2)
{
    return (strcmp(s1, s2) ==0);
}

void fatalError(const char * errorStr)
{
    fprintf(stderr, "%s\n", errorStr);
    exit(1);
}

void fsError(const char * filename)
{
    char * temp;
    if (errno == ENOENT)
        asprintf(&temp, "samblaster: File '%s' does not exist.\n", filename);
    else
        asprintf(&temp, "samblaster: File system error on %s: %d.\n", filename, errno);
    fatalError(temp);
}

///////////////////////////////////////////////////////////////////////////////
// Runtime Statistics
///////////////////////////////////////////////////////////////////////////////

// Stuff needed for timings.
// To turn timing off, set the below to 0.
#define TIMING 1
#if TIMING
// A convenience function for outputing time is seconds in a more useful metric.
void fprintTimeSeconds (FILE * out, double seconds, int precision)
{
    double totalseconds = seconds;
    int hours = seconds/3600.;
    if (hours > 0)
    {
        seconds -= hours * 3600;
        fprintf(out, "%dH", hours);
    }
    int minutes = seconds/60.;
    if (minutes > 0)
    {
        seconds -= minutes * 60;
        fprintf(out, "%dM", minutes);
    }
    if (hours + minutes > 0)
    {
        fprintf(out, "%.0fS", seconds);
        fprintf(out, "(%.*fS)", precision, totalseconds);
    }
    else fprintf(out, "%.*fS", precision, totalseconds);
}

void fprintTimeMicroSeconds (FILE * out, UINT64 microSeconds, int precision)
{
    fprintTimeSeconds(out, ((double)microSeconds/1000000.0), precision);
}

inline UINT64 diffTVs (struct timeval * startTV, struct timeval * endTV)
{
    return (((endTV->tv_sec - startTV->tv_sec) * 1000000) + (endTV->tv_usec - startTV->tv_usec));
}
#include <sys/times.h>
#include <sys/resource.h>
#include <time.h>
#endif // TIMING

///////////////////////////////////////////////////////////////////////////////
// Split Lines
///////////////////////////////////////////////////////////////////////////////

// The structure to store "split" input lines, especially SAM lines.
// They form a singly linked list so that we can form groups of them,
//     and also so that we can keep a freelist of them.

// Reference offsets (pos) can be negative or run off the end of a contig due to clipping.
// Therefore, we will use a padding strategy.
// The space allocated to each contig will be padded by twice the max read length.
// This leaves room for both offset underflow and overflow.
// And all offsets will be shifted higher by the max read length.
// This will eliminate negative offsets and "center" offsets within the offset range for the contig.
typedef INT32 pos_t; 
#define MAX_SEQUENCE_LENGTH 250 // Current illumina paired-end reads are at most 150 + 150
inline int padLength(int length)
{
    return length + (2 * MAX_SEQUENCE_LENGTH);
}
inline int padPos(int pos)
{
    return pos + MAX_SEQUENCE_LENGTH;
}

// We need to pre-define these for the SAM specific fields.
typedef UINT64 sgn_t; // Type for signatures for offsets and lengths.
// And the type itself for the next pointer.
typedef struct splitLine splitLine_t;
splitLine_t * splitLineFreeList = NULL;
struct splitLine
{
    // General fields for any split line.
    splitLine_t * next;
    char * buffer;
    int bufLen;
    size_t maxBufLen;
    char **fields;
    int numFields;
    int maxFields;
    bool split;
    // Special SAM fields that we need to access as other than strings.
    // It this were a class, these would be in a subclass.
    int   flag;
    pos_t pos;
    int   seqNum;
    pos_t binPos;
    int   binNum;
    int   SQO;
    int   EQO;
    int   sclip;
    int   eclip;
    int   rapos;
    int   raLen;
    int   qaLen;
    bool  CIGARprocessed;
    bool  discordant;
    bool  splitter;
    bool  unmappedClipped;
};

// Creator for splitLine
splitLine_t * makeSplitLine()
{
    splitLine_t * line = (splitLine_t *)malloc(sizeof(splitLine_t));
    line->bufLen = 0;
    line->maxBufLen = 1000;
    line->buffer = (char *)malloc(line->maxBufLen);
    line->numFields = 0;
    line->maxFields = 100;
    line->fields = (char **)malloc(line->maxFields * sizeof(char *));
    return line;
}

// Destructor for split line.
void deleteSplitLine(splitLine_t * line)
{
    free(line->buffer);
    free(line->fields);
    free(line);
}

// Descructor for a list of splitLines
void cleanUpSplitLines()
{
    splitLine_t * l = splitLineFreeList;
    while (l != NULL)
    {
        splitLine_t * next = l->next;
        deleteSplitLine(l);
        l = next;
    }
}

// Like descructor for splitLine except don't free memory.
// Instead, put the linked list of objects back on the free list.
void disposeSplitLines(splitLine_t * line)
{
    // First find the last line in the list.
    // Then get rid of them all.
    splitLine_t * last = line;
    for (splitLine_t * l = line->next; l!=NULL; l = l->next) last = l;
    last->next = splitLineFreeList;
    splitLineFreeList = line;
}

// Like constuctor, except take struct off free list if there is one.
splitLine_t * getSplitLine()
{
    splitLine_t * line;
    if (splitLineFreeList ==  NULL)
    {
        line = makeSplitLine();
    }
    else
    {
        line = splitLineFreeList;
        splitLineFreeList = splitLineFreeList->next;
    }
    line->next = NULL;
    line->CIGARprocessed = false;
    // Mark these here so that the code for these doesn't have to stand on its head to do it.
    line->discordant = false;
    line->splitter = false;
    line->unmappedClipped = false;
    return line;
}

// Split the line into fields.
void splitSplitLine(splitLine_t * line, int maxSplits)
{
    line->numFields = 0;
    int fieldStart = 0;
    // replace the newline with a tab so that it works like the rest of the fields.
    line->buffer[line->bufLen-1] = '\t';
    for (int i=0; i<line->bufLen; ++i)
    {
        if (line->buffer[i] == '\t')
        {
            line->fields[line->numFields] = line->buffer + fieldStart;
            line->numFields += 1;
            if (line->numFields == maxSplits) break;
            line->buffer[i] = 0;
            // Get ready for the next iteration.
            fieldStart = i+1;
        }
    }
    // replace the tab at the end of the line with a null char to terminate the final string.
    line->buffer[line->bufLen-1] = 0;
    line->split = true;
}

// Unsplit the fields back into a single string.
// This will mung the strings, so only call this when all processing on the line is done.
void unsplitSplitLine(splitLine_t * line)
{
    // First make sure we are still split.
    if (!line->split) return;
    // First undo the splits.
    // We will undo the splits backwards from the next field to avoid having to calculate strlen each time.
    for (int i=1; i<line->numFields; ++i)
    {
        line->fields[i][-1] = '\t';
    }
    // Now put the newline back in.
    line->buffer[line->bufLen-1] = '\n';
    // Mark as no longer split.
    line->split = false;
}

// Resize the buffer of a splitLine.
// Since the new buffer may not be in the same place, we need to first unsplit, resize, then resplit.
void resizeSplitLine(splitLine_t * line, int newsize)
{
    // First unsplit it.
    unsplitSplitLine(line);
    // Resize the buffer, giving a little extra room.
    line->maxBufLen = newsize + 50;
    line->buffer = (char *)realloc(line->buffer, line->maxBufLen);
    if (line->buffer == NULL)
    {
        fatalError("samblaster: Failed to reallocate to a larger read buffer size.\n");
    }
    // Now resplit the line.
    splitSplitLine(line, line->numFields);
}


// Change a field into a value.
// This will be tough given how we output lines.
// So, we might have to try a few things.
// Start with simply expanding/contracting the string to put in the new value.
void changeFieldSplitLine(splitLine_t * line, int fnum, char * newValue)
{
    // What we do will depend on the lengths of the two strings.
    // So, start by calculaing these only once.
    char * fp = line->fields[fnum];
    int oldLen = strlen(fp);
    int newLen = strlen(newValue);
    // Now see if we need to first change the length of the whole line.
    int move = newLen - oldLen;
    if (move != 0)
    {
        // This should never happen, but to be robust we need to check.
        // It is messy to fix it, as all the field ptrs will now be wrong.
        // For now, punt.
        if ((size_t)(line->bufLen + move) >= line->maxBufLen)
        {
            resizeSplitLine(line, line->bufLen + move);
            fp = line->fields[fnum];
        }
        // Calculate the size of the tail that is still needed.
        int distance = 1 + line->bufLen - (fp - line->buffer) - oldLen;
        // Do the copy.
        memmove(fp+newLen, fp+oldLen, distance);
        // Correct the total length of the buffer.
        line->bufLen += move;
        // We need to correct the other ptrs as well.
        for (int i=fnum+1; i<line->numFields; i++) line->fields[i] += move;
    }
    // Copy in the new value.
    memcpy(fp, newValue, newLen);
}

void addTag(splitLine_t * line, const char * header, const char * val)
{
    int hl = strlen(header);
    int vl = strlen(val);
    // Make sure everything will fit.
    int newlen = line->bufLen + hl + vl;
    if ((size_t)newlen >= line->maxBufLen)
    {
        resizeSplitLine(line, newlen);
    }
    // Copy over the header and the value.
    char * ptr = line->buffer + line->bufLen - 1;
    ptr = (char *)mempcpy(ptr, header, hl);
    ptr = (char *)mempcpy(ptr, val, vl);
    // Add the null terminator for the field, and for the record.
    ptr[0] = 0;
    ptr[1] = 0;
    // Fix the buffer length.
    line->bufLen = newlen;
}


// Read a line from the file and split it.
splitLine_t * readLine(FILE * input)
{
    splitLine_t * sline = getSplitLine();
    sline->bufLen = getline(&sline->buffer, &sline->maxBufLen, input);
    if (sline->bufLen < 1)
    {
        disposeSplitLines(sline);
        return NULL;
    }
    splitSplitLine(sline, 12);
    return sline;
}

inline void outputString(char * str, FILE * output)
{
    // Do the error checking here so we don't have to do it elsewhere.
    if (fputs(str, output) < 0)
    {
        fatalError("samblaster: Unable to write to output file.\n");
    }
}

// Output the line.
inline void writeLine(splitLine_t * line, FILE * output)
{
    unsplitSplitLine(line);
    outputString(line->buffer, output);
}

// Check the first line of a file (e.g. input) for bam signature.
void checkBAMfile(splitLine_t * line)
{
    // If the file is a bam file, we can't rely on fields.
    // So, look at the underlying buffer for the line.

    // First define the signature to look for.
    int values[] = {31, -117, 8, 4, 66, 67,  2};
    int offsets[] = {0,    1, 2, 3, 12, 13, 14};
    int count = 7;

    // Check for empty file or likely a sam file with a header.
    // This is necessary, as the @HD line may not be long enough to check for the BAM signature.
    if (line == NULL) fatalError("samblaster: Input file is empty. Exiting.\n");
    if (line->buffer[0] == '@') return;

    // If a SAM file has no header, an alignment row should easily be long enough.
    if (line->bufLen <= offsets[count-1]) fatalError("samblaster: Input file is empty. Exiting.\n");
    // Check for the BAM signature.
    for (int i=0; i<count; i++)
    {
        if ((int)line->buffer[offsets[i]] != values[i]) return;
    }

    // If we are here, we almost certainly have a bamfile.
    fatalError("samblaster: Input file appears to be in BAM format. SAM input is required. Exiting.\n");
}

///////////////////////////////////////////////////////////////////////////////
// SAM and signature set related structures.
///////////////////////////////////////////////////////////////////////////////

// Define SAM field offsets.
#define QNAME  0
#define FLAG   1
#define RNAME  2
#define POS    3
#define MAPQ   4
#define CIGAR  5
#define RNEXT  6
#define PNEXT  7
#define TLEN   8
#define SEQ    9
#define QUAL  10
#define TAGS  11


// Define SAM flag accessors.
#define MULTI_SEGS     0x1
#define FIRST_SEG      0x40
#define SECOND_SEG     0x80
inline bool checkFlag(splitLine_t * line, int bits) { return ((line->flag & bits) != 0); }

inline void setFlag(splitLine_t * line, int bits) { line->flag |= bits; }

inline bool isPaired(splitLine_t * line) { return checkFlag(line, MULTI_SEGS); }

inline bool isConcordant(splitLine_t * line) { return checkFlag(line, 0x2); }

inline bool isDiscordant(splitLine_t * line) { return !isConcordant(line); }

inline bool isUnmapped(splitLine_t * line) { return checkFlag(line, 0x4); }

inline bool isNextUnmapped(splitLine_t * line) { return checkFlag(line, 0x8); }

inline bool isNextMapped(splitLine_t * line) { return !isNextUnmapped(line); }

inline bool isMapped(splitLine_t * line) { return !isUnmapped(line); }

inline bool isReverseStrand(splitLine_t * line) { return checkFlag(line, 0x10); }

inline bool isForwardStrand(splitLine_t * line) { return !isReverseStrand(line); }

inline bool isFirstRead(splitLine_t * line) { return checkFlag(line, FIRST_SEG); }

inline bool isSecondRead(splitLine_t * line) { return checkFlag(line, SECOND_SEG); }

// These determine alignment type.
// Things may get more complicated than this once we have alternate contigs such as in build 38 of human genome.
inline bool isPrimaryAlignment(splitLine_t * line)
{ return !(checkFlag(line, 0x100) || checkFlag(line, 0x800)); }

// We have to hande secondard and complementary alignments differently depending on compatMode.
// So, we store which bits are being included in each.

int complementaryBits = 0x800;
inline bool isComplementaryAlignment(splitLine_t * line)
{ return checkFlag(line, complementaryBits); }

int secondaryBits = 0x100;
inline bool isSecondaryAlignment(splitLine_t * line)
{ return checkFlag(line, secondaryBits); }

inline bool isDuplicate(splitLine_t * line) { return checkFlag(line, 0x400); }

inline void setDuplicate(splitLine_t * line) { setFlag(line, 0x400); }

typedef hashTable_t sigSet_t;

inline int str2int (char * str)
{
    return strtol(str, NULL, 0);
}

// Need to change this if pos is unsigned.
inline pos_t str2pos (char * str)
{
    return strtol(str, NULL, 0);
}

// Temp buffer to use to form new flag field when marking dups.
char tempBuf[10];
inline void markDup(splitLine_t * line)
{
    setDuplicate(line);
    sprintf(tempBuf, "%d", line->flag);
    changeFieldSplitLine(line, FLAG, tempBuf);
}

// Special version of write line that appends an id number to the output.
// Used to output splitters.
void writeSAMlineWithIdNum(splitLine_t * line, FILE * output)
{
    // Unsplit the line.
    unsplitSplitLine(line);
    // Split it two ways to isolate the id field.
    splitSplitLine(line, 2);
    outputString(line->fields[0], output);
    if (isPaired(line)) fprintf(output, "_%d\t", isFirstRead(line) ? 1 : 2);
    else                fprintf(output, "\t");
    outputString(line->fields[1], output);
    fprintf(output, "\n");
}

///////////////////////////////////////////////////////////////////////////////
// Sequence Map
///////////////////////////////////////////////////////////////////////////////

// We use a map instead of a hash map for sequence names.
// This is because the default hash function on char * hashes the ptr values.
// So, we would need to define our own hash on char * to get things to work properly.
// Not worth it for a structure holding so few members.

// Function needed to get char * map to work.
struct less_str
{
   bool operator()(char const *a, char const *b) const
   {
      return strcmp(a, b) < 0;
   }
};

// This stores the map between sequence names and sequence numbers.
typedef std::map<const char *, int, less_str> seqMap_t;

inline void addSeq(seqMap_t * seqs, char * item, int val)
{
    (*seqs)[item] = val;
}

///////////////////////////////////////////////////////////////////////////////
// Struct for processing state
///////////////////////////////////////////////////////////////////////////////
struct state_struct
{
    char *         inputFileName;
    FILE *         inputFile;
    char *         outputFileName;
    FILE *         outputFile;
    FILE *         discordantFile;
    char *         discordantFileName;
    FILE *         splitterFile;
    char *         splitterFileName;
    FILE *         unmappedClippedFile;
    char *         unmappedClippedFileName;
    sigSet_t *     sigs;
    seqMap_t       seqs;
    UINT32 *       seqLens;
    UINT64 *       seqOffs;
    splitLine_t ** splitterArray;
    int            splitterArrayMaxSize;
    UINT32         sigArraySize;
    int            binCount;
    int            minNonOverlap;
    int            maxSplitCount;
    int            minIndelSize;
    int            maxUnmappedBases;
    int            minClip;
    int            unmappedFastq;
    bool           acceptDups;
    bool           excludeDups;
    bool           removeDups;
    bool           addMateTags;
    bool           compatMode;
    bool           ignoreUnmated;
    bool           quiet;
};
typedef struct state_struct state_t;

state_t * makeState ()
{
    state_t * s = new state_t();
    s->inputFile = stdin;
    s->inputFileName = (char *)"stdin";
    s->outputFile = stdout;
    s->outputFileName = (char *)"stdout";
    s->discordantFile = NULL;
    s->discordantFileName = (char *)"";
    s->splitterFile = NULL;
    s->splitterFileName = (char *)"";
    s->unmappedClippedFile = NULL;
    s->unmappedClippedFileName = (char *)"";
    s->sigs = NULL;
    s->minNonOverlap = 20;
    s->maxSplitCount = 2;
    s->minIndelSize = 50;
    s->maxUnmappedBases = 50;
    s->minClip = 20;
    s->acceptDups = false;
    s->excludeDups = false;
    s->removeDups = false;
    s->addMateTags = false;
    s->compatMode = false;
    s->ignoreUnmated = false;
    s->quiet = false;
    // Start this as -1 to indicate we don't know yet.
    // Once we are outputting our first line, we will decide.
    s->unmappedFastq = -1;
    // Used as a temporary location for ptrs to splitter for sort routine.
    s->splitterArrayMaxSize = 1000;
    s->splitterArray = (splitLine_t **)(malloc(s->splitterArrayMaxSize * sizeof(splitLine_t *)));
    return s;
}

void deleteState(state_t * s)
{
    free(s->splitterArray);
    if (s->sigs != NULL)
    {
        // delete[] s->sigs;
        for (UINT32 i=0; i<s->sigArraySize; i++) deleteHashTable(&(s->sigs[i]));
        free (s->sigs);
    }
    for (seqMap_t::iterator iter = s->seqs.begin(); iter != s->seqs.end(); ++iter)
    {
        free((char *)(iter->first));
    }
    if (s->seqLens != NULL) free(s->seqLens);
    if (s->seqOffs != NULL) free(s->seqOffs);
    delete s;
}


///////////////////////////////////////////////////////////////////////////////
// Signatures
///////////////////////////////////////////////////////////////////////////////
// We now calculate signatures as offsets into a super contig that includes the entire genome.
// And partition it into equaly sized bins.
// This performs better on genomes with large numbers of small contigs,
//  without performance degradation on more standard genomes.
// Thanks to https://github.com/carsonhh for the suggestion.
/////////////////////////////////////////////////////////////////////////////
inline sgn_t calcSig(splitLine_t * first, splitLine_t * second)
{
    // Total nonsense to get the compiler to actually work.
    UINT64 t1 = first->binPos;
    UINT64 t2 = t1 << 32;
    UINT64 final = t2 | second->binPos;
    return (sgn_t)final;
}

inline UINT32 calcSigArrOff(splitLine_t * first, splitLine_t * second, int binCount)
{
    UINT32 s1 = (first->binNum * 2) + (isReverseStrand(first) ? 1 : 0);
    UINT32 s2 = (second->binNum * 2) + (isReverseStrand(second) ? 1 : 0);
    UINT32 retval = (s1 * binCount * 2) + s2;
#ifdef DEBUG
    fprintf(stderr, "1st %d %d -> %d 2nd %d %d -> %d count %d result %d read: %s\n",
            first->binNum, isReverseStrand(first), s1, second->binNum, isReverseStrand(second), s2, binCount, retval, first->fields[QNAME]);
#endif
    return retval;
}

///////////////////////////////////////////////////////////////////////////////
// Sequences
///////////////////////////////////////////////////////////////////////////////

inline int getSeqNum(splitLine_t * line, int field, state_t * state) __attribute__((always_inline));
inline int getSeqNum(splitLine_t * line, int field, state_t * state)
{
#ifdef DEBUG
    seqMap_t::iterator findret = state->seqs.find(line->fields[field]);
    if (findret == state->seqs.end())
    {
        char * temp;
        asprintf(&temp, "Unable to find seq %s for readid %s in sequence map.\n", line->fields[field], line->fields[QNAME]);
        fatalError(temp);
    }
    return findret->second;
#else
    return state->seqs.find(line->fields[field])->second;
#endif
}

///////////////////////////////////////////////////////////////////////////////
// Helpers to process CIGAR strings
///////////////////////////////////////////////////////////////////////////////

// This will parse a base 10 int, and change ptr to one char beyond the end of the number.
inline int parseNextInt(char **ptr)
{
    int num = 0;
    for (char curChar = (*ptr)[0]; curChar != 0; curChar = (++(*ptr))[0])
    {
        int digit = curChar - '0';
        if (digit >= 0 && digit <= 9) num = num*10 + digit;
        else break;
    }
    return num;
}

// This will the current char, and move the ptr ahead by one.
inline char parseNextOpCode(char **ptr)
{
    return ((*ptr)++)[0];
}

// This just test for end of string.
inline bool moreCigarOps(char *ptr)
{
    return (ptr[0] != 0);
}

void calcOffsets(splitLine_t * line)
{
    if (line->CIGARprocessed) return;
    char * cigar = line->fields[CIGAR];
    line->raLen = 0;
    line->qaLen = 0;
    line->sclip = 0;
    line->eclip = 0;
    bool first = true;
    while (moreCigarOps(cigar))
    {
        int opLen = parseNextInt(&cigar);
        char opCode = parseNextOpCode(&cigar);
        if      (opCode == 'M' || opCode == '=' || opCode == 'X')
        {
            line->raLen += opLen;
            line->qaLen += opLen;
            first = false;
        }
        else if (opCode == 'S' || opCode == 'H')
        {
            if (first) line->sclip += opLen;
            else       line->eclip += opLen;
        }
        else if (opCode == 'D' || opCode == 'N')
        {
            line->raLen += opLen;
        }
        else if (opCode == 'I')
        {
            line->qaLen += opLen;
        }
        else
        {
            fprintf(stderr, "Unknown opcode '%c' in CIGAR string: '%s'\n", opCode, line->fields[CIGAR]);
        }
    }
    line->rapos = str2pos(line->fields[POS]);
    if (isForwardStrand(line))
    {
        line->pos = line->rapos - line->sclip;
        line->SQO = line->sclip;
        line->EQO = line->sclip + line->qaLen - 1;
    }
    else
    {
        line->pos = line->rapos + line->raLen + line->eclip - 1;
        line->SQO = line->eclip;
        line->EQO = line->eclip + line->qaLen - 1;
    }
    // Need to pad the pos in case it is negative
    line->pos = padPos(line->pos);
    // Let's not calculate these again for this line.
    line->CIGARprocessed = true;
}

inline int getStartDiag(splitLine_t * line)
{
    // SRO - SQO (not strand normalized)
    // Simplify the following.
    // return (str2pos(line->fields[POS])) - line->sclip;
    return line->rapos - line->sclip;
}

inline int getEndDiag(splitLine_t * line)
{
    // ERO - EQO (not strand normalized)
    // Simplify the following
    // return (line->rapos + line->raLen - 1) - (line->sclip + line->qaLen - 1)
    return (line->rapos + line->raLen) - (line->sclip + line->qaLen);
}


///////////////////////////////////////////////////////////////////////////////
// Process SAM Blocks
///////////////////////////////////////////////////////////////////////////////
#define BIN_SHIFT 27 //bin window is 27 bits wide
#define BIN_MASK ((1 << 27)-1) //bin window is 27 bits wide

// This is apparently no longer called.
void outputSAMBlock(splitLine_t * block, FILE * output)
{
    for (splitLine_t * line = block; line != NULL; line = line->next)
    {
        writeLine(line, output);
    }
    disposeSplitLines(block);
}

inline bool needSwap(splitLine_t * first, splitLine_t * second)
{
    // Sort first by ref offset.
    if (first->pos > second->pos) return true;
    if (first->pos < second->pos) return false;
    // Now by seq number.
    if (first->seqNum > second->seqNum) return true;
    if (first->seqNum < second->seqNum) return false;
    // Now by strand.
    // If they are both the same strand, it makes no difference which on is first.
    if (isReverseStrand(first) == isReverseStrand(second)) return false;
    if (isReverseStrand(first) && isForwardStrand(second)) return true;
    return false;
}

inline void swapPtrs(splitLine_t ** first, splitLine_t ** second)
{
    splitLine_t * temp = *first;
    *first = *second;
    *second = temp;
}

void brokenBlock(splitLine_t *block, int count)
{
    char * temp;
    asprintf(&temp, "samblaster: Can't find first and/or second of pair in sam block of length %d for id: %s\n%s%s:%s\n%s",
             count, block->fields[QNAME], "samblaster:    At location: ", block->fields[RNAME], block->fields[POS],
             "samblaster:    Are you sure the input is sorted by read ids?");
    fatalError(temp);
}

// Some fields for statistics.
UINT64 idCount = 0;
UINT64 dupCount = 0;
UINT64 discCount = 0;
UINT64 splitCount = 0;
UINT64 unmapClipCount = 0;
UINT64 unmatedCount = 0;

// This is the main workhorse that determines if lines are dups or not.
template <bool excludeSecondaries>
int fillSplitterArray(splitLine_t * block, state_t * state, int mask, bool flagValue);
void markDupsDiscordants(splitLine_t * block, state_t * state)
{
    splitLine_t * first = NULL;
    splitLine_t * second = NULL;
    int count = 0;
    for (splitLine_t * line = block; line != NULL; line = line->next)
    {
        count += 1;
        // Do this conversion once and store the result.
        line->flag = str2int(line->fields[FLAG]);
        // We make our duplicate decisions based solely on primary alignments.
        if (!isPrimaryAlignment(line)) continue;
        // Allow unpaired reads to go through (as the second so that signature is correct).
        // According to the SAM spec, this must be checked first.
        if (!isPaired(line)) second = line;
        // Figure out if this is the first half or second half of a pair.
        else if (isFirstRead(line))  first  = line;
        else if (isSecondRead(line)) second = line;
    }

    // Figure out what type of "pair" we have.
    bool orphan = false;
    bool dummyFirst = false;
    // First get rid of the useless case of having no first AND no second.
    if (first == NULL && second == NULL) goto outOfHere;
    // Now see if we have orphan with the unmapped read missing.
    if (first == NULL || second == NULL)
    {
        // Get the NULL one in the first slot.
        if (second == NULL) swapPtrs(&first, &second);
        // If the only read says its paired, and it is unmapped or its mate is mapped, something is wrong.
        if (isPaired(second) && (isUnmapped(second) || isNextMapped(second))) goto outOfHere;
        // If the only read we have is unmapped, then it can't be a dup.
        if (isUnmapped(second)) return;
        // Now MAKE a dummy record for the first read, but don't put it into the block.
        // That way we won't have to worry about it getting output, or when processing splitters etc.
        // But, we will have to remember to dispose of it.
        first = getSplitLine();
        // Set the flag field to what it would have been.
        // What if this "pair" is a singleton read?
        first->flag = isFirstRead(second) ? 0x85 : 0x45;
        orphan = true;
        dummyFirst = true;
    }
    else
    {
        // Handle the addition of MC and MQ tags if requested.
        if (state->addMateTags)
        {
            int mask = (FIRST_SEG | SECOND_SEG);
            // Process the first of the pair.
            // Get the list of reads that match the second of the pair.
            if (isMapped(first))
            {
                int count = fillSplitterArray<false>(block, state, second->flag & mask, true);
                for (int i=0; i<count; ++i)
                {
                    splitLine_t * line = state->splitterArray[i];
                    addTag(line, "	MC:Z:", first->fields[CIGAR]);
                    addTag(line, "	MQ:i:", first->fields[MAPQ]);
                }
            }
            // Process the second of the pair.
            // Get the list of reads that match the first of the pair.
            if (isMapped(second))
            {
                count = fillSplitterArray<false>(block, state, first->flag & mask, true);
                for (int i=0; i<count; ++i)
                {
                    splitLine_t * line = state->splitterArray[i];
                    addTag(line, "	MC:Z:", second->fields[CIGAR]);
                    addTag(line, "	MQ:i:", second->fields[MAPQ]);
                }
            }
        }

        // Never mark pairs as dups if both sides are unmapped.
        if (isUnmapped(first) && isUnmapped(second)) return;

        // We need to properly handle orphans to get the correct reference offsets and sequence numbers.
        orphan = (isUnmapped(first) || isUnmapped(second));
        // Orphan that needs to be swapped.
        // We need the unmapped one in the first slot so that they won't all collide in the hash table.
        if (isMapped(first) && isUnmapped(second))
        {
            swapPtrs(&first, &second);
        }
    }

    // Now look for duplicates.
    if (!state->acceptDups)
    {
        // Calculate and store the second position and sequence name.
        calcOffsets(second);
        second->seqNum = getSeqNum(second, RNAME, state);
        UINT64 seqOff = state->seqOffs[second->seqNum]; //genome relative position
        second->binNum = (seqOff + second->pos) >> BIN_SHIFT;
        second->binPos = (seqOff + second->pos) &  BIN_MASK;

        if (orphan)
        {
            // We have an orphan, so we just zero out the pos and seqnum
            first->pos = 0;
            first->seqNum = 0;
            first->binNum = 0;
            first->binPos = 0;
        }
        else
        {
            // Not an orphan, so handle first on its own.
            calcOffsets(first);
            first->seqNum = getSeqNum(first, RNAME, state);
            seqOff = state->seqOffs[first->seqNum]; //genome relative position
            first->binNum = (seqOff + first->pos) >> BIN_SHIFT;
            first->binPos = (seqOff + first->pos) &  BIN_MASK;
        }

        // The fact of which alignment is first or second in the template is not relevant for determining dups.
        // Therefore, we normalize the pairs based on their characteristics.
        // We have already swapped orphans.
        // Otherwise, sort by pos, and if equal, sort by sequence num, then by strand.
        if (!orphan && needSwap(first, second)) swapPtrs(&first, &second);

        // Now find the signature of the pair.
        sgn_t sig = calcSig(first, second);
        // Calculate the offset into the signatures array.
        UINT32 off = calcSigArrOff(first, second, state->binCount);
        // Attempt insert into the sigs structure.
        // The return value will tell us if it was already there.
        bool insert = hashTableInsert(&(state->sigs[off]), sig);
        // Check if the insertion actually happened.
        if (!insert)
        {
            dupCount += 1;
            // We always mark all or none of a block as dup.
            for (splitLine_t * line = block; line != NULL; line = line->next)
            {
                markDup(line);
            }
        }
    }

    // If we have a dummy first, we can't have a discordant pair.
    if (dummyFirst)
    {
        disposeSplitLines(first);
        return;
    }

    // The first and second help us mark the discordants.
    // Both sides mapped, but pair not properly aligned.
    if (!orphan && isDiscordant(first))
    {
        first->discordant = true;
        second->discordant = true;
    }
    return;

outOfHere:
    if (state->ignoreUnmated) {unmatedCount += 1; return;}
    else                       brokenBlock(block, count);
}

// Sort ascending in SQO.
int compQOs(const void * p1, const void * p2)
{
    splitLine_t * l1 = (*(splitLine_t **)p1);
    splitLine_t * l2 = (*(splitLine_t **)p2);
    return (l1->SQO - l2->SQO);
}

template <bool excludeSecondaries>
int fillSplitterArray(splitLine_t * block, state_t * state, int mask, bool flagValue)
{
    // Count the secondaries we have for this read (if any), and store their ptrs into an array.
    int count = 0;
    for (splitLine_t * line = block; line != NULL; line = line->next)
    {
        // For all the ones that are the current read of interest....
        // Check if they are a primary or complementary alignment.
        if (checkFlag(line, mask) == flagValue && !(excludeSecondaries && isSecondaryAlignment(line)))
        {
            // Add the ptr to this line to the sort array.
            // If it won't fit, double the array size.
            if (count >= state->splitterArrayMaxSize)
            {
                state->splitterArrayMaxSize *= 2;
                state->splitterArray = (splitLine_t **)(realloc(state->splitterArray,
                                                                state->splitterArrayMaxSize * sizeof(splitLine_t *)));
            }
            state->splitterArray[count] = line;
            count += 1;
        }
    }
    return count;
}

void markSplitterUnmappedClipped(splitLine_t * block, state_t * state, int mask, bool flagValue)
{
    // Count the secondaries we have for this read (if any), and store their ptrs into an array.
    int count = fillSplitterArray<true>(block, state, mask, flagValue);

    // We have the lines of interest in an array.
    // Decide what to do next based on the number of reads.
    if (count == 0) return;
    if (count == 1)
    {
        if (state->unmappedClippedFile == NULL) return;
        // Process unmapped or clipped.
        splitLine_t * line = state->splitterArray[0];
        // Unmapped or clipped alignments should be primary.
        if (!isPrimaryAlignment(line)) return;
        if (isUnmapped(line))
        {
            line->unmappedClipped = true;
        }
        else
        {
            // Process the CIGAR string.
            // As this is expensive, we delay as long as possible.
            calcOffsets(line);
            if (line->sclip >= state->minClip || line->eclip >= state->minClip)
            {
                line->unmappedClipped = true;
            }
        }
        return;
    }

    // See if we need to process for splitters.
    if (state->splitterFile == NULL || count > state->maxSplitCount) return;

    // Calculate the query positions (for sorting) and do other preprocessing.
    for (int i=0; i<count; i++)
    {
        splitLine_t * line = state->splitterArray[i];
        // Make sure the primary is mapped!
        if (isPrimaryAlignment(line) && isUnmapped(line)) return;
        calcOffsets(line);
    }

    // We need to sort it by strand normalized query offsets.
    qsort(state->splitterArray, count, sizeof(splitLine_t *), compQOs);

    // Now check for pairs that match the desired parameters.
    splitLine_t * left = state->splitterArray[0];
    splitLine_t * right;
    for (int i=1; i<count; i++)
    {
        right = state->splitterArray[i];

        // First check for minNonOverlap.
        // We don't allow negative overlap, as that will lead to wrong non-overlap calculation.
        int overlap = std::max(1 + std::min(left->EQO, right->EQO) - std::max(left->SQO, right->SQO), 0);
        int alen1 = 1 + left->EQO - left->SQO;
        int alen2 = 1 + right->EQO - right->SQO;
        int mno = std::min(alen1-overlap, alen2-overlap);
        if (mno < state->minNonOverlap) goto nextIter;

        // Now check for the deserts and diagonal difference.
        // If they are on different chroms or strands, they pass without the other checks.
        // Since we only care if the sequences are the same, we don't need seqNums.
        // Instead just compare the strings!
        if (streq(left->fields[RNAME], right->fields[RNAME]) && (isReverseStrand(left) == isReverseStrand(right)))
        {
            // The start and end diags might be different if there is an indel in the alignments.
            // So, we use the end on the left, and the start on the right.
            // This will give us the net diag difference between them.
            int leftDiag, rightDiag, insSize;
            if (isReverseStrand(left))
            {
                leftDiag = getStartDiag(left);
                rightDiag = getEndDiag(right);
                insSize = rightDiag - leftDiag;
            }
            else
            {
                leftDiag = getEndDiag(left);
                rightDiag = getStartDiag(right);
                insSize = leftDiag - rightDiag;
            }
            int desert = right->SQO - left->EQO - 1;
            // The absolute value will handle both inserts and deletes.
            // So check indel is big enough, and that there are not too many unmapped bases.
            // Subtract the inadvertant desert gap of inserts.
            if ((abs(insSize) < state->minIndelSize)
                || ((desert > 0) && ((desert - (int)std::max(0, insSize)) > state->maxUnmappedBases)))
                goto nextIter;
        }
        // We made it through the gamet.
        // Mark this pair for output.
        left->splitter = true;
        right->splitter = true;

    nextIter:
        // Get ready for the next iteration.
        left = right;
    }
}

void writeUnmappedClipped(splitLine_t * line, state_t * state)
{
    // Check if we are outputting fasta or fastq.
    if (state->unmappedFastq == -1)
        state->unmappedFastq = (streq(line->fields[QUAL], "*") ? 0 : 1);

    // Print the first line.
    char firstChar = (state->unmappedFastq) ? '@' : '>';
    if (isPaired(line)) fprintf(state->unmappedClippedFile, "%c%s_%d\n", firstChar, line->fields[QNAME], (isFirstRead(line) ? 1 : 2));
    else                fprintf(state->unmappedClippedFile, "%c%s\n", firstChar, line->fields[QNAME]);

    if (state->unmappedFastq) fprintf(state->unmappedClippedFile, "%s\n+\n%s\n", line->fields[SEQ], line->fields[QUAL]);
    else                       fprintf(state->unmappedClippedFile, "%s\n", line->fields[SEQ]);
}

void processSAMBlock(splitLine_t * block, state_t * state)
{
    idCount += 1;
    // First mark dups and find the discordants.
    // These share a lot of looking at flag bits, so make sense to put together.
    markDupsDiscordants(block, state);
    // Look for splitters.
    // Since this is expensive, don't do it unless the user asked for them.
    if (state->splitterFile != NULL || state->unmappedClippedFile != NULL)
    {
        // Check the first read for splitter.
        markSplitterUnmappedClipped(block, state, FIRST_SEG, true);
        // Check the second read for splitter.
        markSplitterUnmappedClipped(block, state, SECOND_SEG, true);
        // Check for a singleton read
        markSplitterUnmappedClipped(block, state, MULTI_SEGS, false);
    }

    // Now do the output.
    for (splitLine_t * line = block; line != NULL; line = line->next)
    {
        // Do the unmapped file first, as it is not sam, and doesn't sew the line back together.
        if (state->unmappedClippedFile != NULL && line->unmappedClipped && !(state->excludeDups && isDuplicate(line)))
        {
            writeUnmappedClipped(line, state);
            unmapClipCount += 1;
        }

        // Write to the output file.
        if (!(state->removeDups && isDuplicate(line)))
        {
            writeLine(line, state->outputFile);
        }

        // Write to discordant file if appropriate.
        if (state->discordantFile != NULL && line->discordant && !(state->excludeDups && isDuplicate(line)))
        {
            writeLine(line, state->discordantFile);
            discCount += 1;
        }
        // Write to splitter file if appropriate.
        if (state->splitterFile != NULL && line->splitter && !(state->excludeDups && isDuplicate(line)))
        {
            writeSAMlineWithIdNum(line, state->splitterFile);
            splitCount += 1;
        }
    }

    disposeSplitLines(block);
}


///////////////////////////////////////////////////////////////////////////////
// Main Routine with helpers.
///////////////////////////////////////////////////////////////////////////////
void printPGsamLine(FILE * f, state_t * s)
{
    if (f == NULL) return;
    fprintf(f, "@PG\tID:SAMBLASTER\tVN:0.1.%d\tCL:samblaster -i %s -o %s", BUILDNUM, s->inputFileName, s->outputFileName);
    if (s->compatMode) fprintf(f, " -M");
    if (s->acceptDups) fprintf(f, " --acceptDupMarks");
    if (s->removeDups) fprintf(f, " --removeDups");
    else if (s->excludeDups && (s->discordantFile != NULL || s->splitterFile != NULL || s->unmappedClippedFile != NULL)) fprintf(f, " --excludeDups");
    if (s->addMateTags) fprintf(f, " --addMateTags");
    if (s->ignoreUnmated) fprintf(f, " --ignoreUnmated");
    if (s->discordantFile != NULL) fprintf(f, " -d %s", s->discordantFileName);
    if (s->splitterFile != NULL)
        fprintf(f, " -s %s --maxSplitCount %d --maxUnmappedBases %d --minIndelSize %d --minNonOverlap %d",
                s->splitterFileName, s->maxSplitCount, s->maxUnmappedBases, s->minIndelSize, s->minNonOverlap);
    if (s->unmappedClippedFile != NULL)
        fprintf(f, " -u %s --minClipSize %d", s->unmappedClippedFileName, s->minClip);
    fprintf(f, "\n");
}

void printVersionString()
{
    fprintf(stderr, "samblaster: Version 0.1.%d\n", BUILDNUM);
}

void printUsageString()
{
    const char* useString =
        "Author: Greg Faust (gf4ea@virginia.edu)\n"
        "Tool to mark duplicates and optionally output split reads and/or discordant pairs.\n"
        "Input sam file must contain paired end data, contain sequence header and be sorted by read ids.\n"
        "Output will be all alignments in the same order as input, with duplicates marked with FLAG 0x400.\n\n"

        "Usage:\n"
        "For use as a post process on an aligner (eg. bwa mem):\n"
        "     bwa mem <idxbase> samp.r1.fq samp.r2.fq | samblaster [-e] [-d samp.disc.sam] [-s samp.split.sam] | samtools view -Sb - > samp.out.bam\n"
        "     bwa mem -M <idxbase> samp.r1.fq samp.r2.fq | samblaster -M [-e] [-d samp.disc.sam] [-s samp.split.sam] | samtools view -Sb - > samp.out.bam\n"
        "For use with a pre-existing bam file to pull split, discordant and/or unmapped reads:\n"
        "     samtools view -h samp.bam | samblaster [-a] [-e] [-d samp.disc.sam] [-s samp.split.sam] [-u samp.umc.fasta] -o /dev/null\n\n"

        "Input/Output Options:\n"
        "-i --input           FILE Input sam file [stdin].\n"
        "-o --output          FILE Output sam file for all input alignments [stdout].\n"
        "-d --discordantFile  FILE Output discordant read pairs to this file. [no discordant file output]\n"
        "-s --splitterFile    FILE Output split reads to this file abiding by paramaters below. [no splitter file output]\n"
        "-u --unmappedFile    FILE Output unmapped/clipped reads as FASTQ to this file abiding by parameters below. [no unmapped file output].\n"
        "                          Requires soft clipping in input file.  Will output FASTQ if QUAL information available, otherwise FASTA.\n\n"

        "Other Options:\n"
        "-a --acceptDupMarks       Accept duplicate marks already in input file instead of looking for duplicates in the input.\n"
        "-e --excludeDups          Exclude reads marked as duplicates from discordant, splitter, and/or unmapped file.\n"
        "-r --removeDups           Remove duplicates reads from all output files. (Implies --excludeDups).\n"
        "   --addMateTags          Add MC and MQ tags to all output paired-end SAM lines.\n"
        "   --ignoreUnmated        Suppress abort on unmated alignments. Use only when sure input is read-id grouped and alignments have been filtered.\n"
        "-M                        Run in compatibility mode; both 0x100 and 0x800 are considered chimeric. Similar to BWA MEM -M option.\n"
        "   --maxSplitCount    INT Maximum number of split alignments for a read to be included in splitter file. [2]\n"
        "   --maxUnmappedBases INT Maximum number of un-aligned bases between two alignments to be included in splitter file. [50]\n"
        "   --minIndelSize     INT Minimum structural variant feature size for split alignments to be included in splitter file. [50]\n"
        "   --minNonOverlap    INT Minimum non-overlaping base pairs between two alignments for a read to be included in splitter file. [20]\n"
        "   --minClipSize      INT Minumum number of bases a mapped read must be clipped to be included in unmapped file. [20]\n"
        "-q --quiet                Output fewer statistics.\n";

        printVersionString();
        fprintf(stderr, "%s", useString);
}

void printUsageStringAbort()
{
    printUsageString();
    exit(1);
}

int main (int argc, char *argv[])
{
    // Set up timing.
#if TIMING
    time_t startTime;
    startTime = time(NULL);
    struct rusage usagebuf;
    getrusage(RUSAGE_SELF, &usagebuf);
    struct timeval startRUTime;
    startRUTime = usagebuf.ru_utime;
#endif

    // Parse input parameters.
    state_t * state = makeState();
    for (int argi = 1; argi< argc; argi++)
    {
        // First the general mode and I/O parameters.
        if (streq(argv[argi], "-h") || streq(argv[argi], "--help"))
        {
            printUsageString();
            return 0;
        }
        else if (streq(argv[argi], "--version"))
        {
            printVersionString();
            return 0;
        }
        else if (streq(argv[argi], "-a") || streq(argv[argi],"--acceptDupMarks"))
        {
            state->acceptDups = true;
        }
        else if (streq(argv[argi], "-d") || streq(argv[argi],"--discordantFile"))
        {
            argi++;
            state->discordantFileName = argv[argi];
        }
        else if (streq(argv[argi], "-e") || streq(argv[argi],"--excludeDups"))
        {
            state->excludeDups = true;
        }
        else if (streq(argv[argi], "-i") || streq(argv[argi],"--input"))
        {
            argi++;
            state->inputFileName = argv[argi];
            if (streq(state->inputFileName, "-")) state->inputFileName = (char *)"stdin";
        }
        else if (streq(argv[argi],"--addMateTags"))
        {
            state->addMateTags = true;
        }
        else if (streq(argv[argi],"--ignoreUnmated"))
        {
            state->ignoreUnmated = true;
        }
        else if (streq(argv[argi],"-M"))
        {
            state->compatMode = true;
            // In compatibility mode, both 0x100 and 0x800 are considered chimeric and can be used as splitters.
            // None will be secondary.
            complementaryBits = (0x100 | 0x800);
            secondaryBits = 0;
        }
        else if (streq(argv[argi],"--maxSplitCount"))
        {
            argi++;
            state->maxSplitCount = str2int(argv[argi]);
        }
        else if (streq(argv[argi],"--maxUnmappedBases"))
        {
            argi++;
            state->maxUnmappedBases = str2int(argv[argi]);
        }
        else if (streq(argv[argi],"--minClipSize"))
        {
            argi++;
            state->minClip = str2int(argv[argi]);
        }
        else if (streq(argv[argi],"--minIndelSize"))
        {
            argi++;
            state->minIndelSize = str2int(argv[argi]);
        }
        else if (streq(argv[argi],"--minNonOverlap"))
        {
            argi++;
            state->minNonOverlap = str2int(argv[argi]);
        }
        else if (streq(argv[argi], "-o") || streq(argv[argi],"--output"))
        {
            argi++;
            state->outputFileName = argv[argi];
        }
        else if (streq(argv[argi], "-r") || streq(argv[argi],"--removeDups"))
        {
            state->removeDups = true;
        }
        else if (streq(argv[argi], "-q") || streq(argv[argi], "--quiet"))
        {
            state->quiet = true;
        }
        else if (streq(argv[argi], "-s") || streq(argv[argi],"--splitterFile"))
        {
            argi++;
            state->splitterFileName = argv[argi];
        }
        else if (streq(argv[argi], "-u") || streq(argv[argi],"--unmappedFile"))
        {
            argi++;
            state->unmappedClippedFileName = argv[argi];
        }
        else
        {
            fprintf(stderr, "samblaster: Unrecognized option: %s\n", argv[argi]);
            printUsageString();
            return 1;
        }
    }

    // Make removeDups imply excludeDups
    if (state->removeDups) state->excludeDups = true;

    // Output the version number.
    printVersionString();

    // Error check and open files.
    // Start with the input.
    if (state->splitterFile != NULL && state->minNonOverlap < 1)
    {
        char * temp;
        asprintf(&temp, "samblaster: Invalid minimum non overlap parameter given: %d\n", state->minNonOverlap);
        fatalError(temp);
    }
    if (state->splitterFile != NULL && state->maxSplitCount < 2)
    {
        char * temp;
        asprintf(&temp, "samblaster: Invalid maximum split count parameter given: %d\n", state->maxSplitCount);
        fatalError(temp);
    }
    if (streq(state->inputFileName, "stdin"))
    {
        fprintf(stderr, "samblaster: Inputting from stdin\n");
    }
    else
    {
        if (!state->quiet) fprintf(stderr, "samblaster: Opening %s for read.\n", state->inputFileName);
        state->inputFile = fopen(state->inputFileName, "r");
        if (state->inputFile == NULL) fsError(state->inputFileName);
    }
    if (streq(state->outputFileName, "stdout"))
    {
        if (!state->quiet) fprintf(stderr, "samblaster: Outputting to stdout\n");
    }
    else
    {
        if (!state->quiet) fprintf(stderr, "samblaster: Opening %s for write.\n", state->outputFileName);
        state->outputFile = fopen(state->outputFileName, "w");
        if (state->outputFile == NULL) fsError(state->outputFileName);
    }
    if (!streq(state->discordantFileName, ""))
    {
        if (!state->quiet) fprintf(stderr, "samblaster: Opening %s for write.\n", state->discordantFileName);
        state->discordantFile = fopen(state->discordantFileName, "w");
        if (state->discordantFile == NULL) fsError(state->discordantFileName);
    }
    if (!streq(state->splitterFileName, ""))
    {
        if (!state->quiet) fprintf(stderr, "samblaster: Opening %s for write.\n", state->splitterFileName);
        state->splitterFile = fopen(state->splitterFileName, "w");
        if (state->splitterFile == NULL) fsError(state->splitterFileName);
    }
    if (!streq(state->unmappedClippedFileName, ""))
    {
        if (!state->quiet) fprintf(stderr, "samblaster: Opening %s for write.\n", state->unmappedClippedFileName);
        state->unmappedClippedFile = fopen(state->unmappedClippedFileName, "w");
        if (state->unmappedClippedFile == NULL) fsError(state->unmappedClippedFileName);
    }

    // Read in the SAM header and create the seqs structure.
    state->seqLens = (UINT32*)calloc(1, sizeof(UINT32)); //initialize to 0
    state->seqOffs = (UINT64*)calloc(1, sizeof(UINT64)); //initialize to 0
    state->seqs[strdup("*")] = 0;
    state->seqLens[0] = padLength(0);
    state->seqOffs[0] = 0;
    int count = 1;
    UINT64 totalLen = 0;
    splitLine_t * line;
    // Read the first line to prime the loop, and also to allow checking for malformed input.
    line = readLine(state->inputFile);
    checkBAMfile(line);
    while (line != NULL && line->fields[0][0] == '@')
    {
        // Process input line to see if it defines a sequence.
        if (streq(line->fields[0], "@SQ"))
        {
            char * seqID = NULL;
            int seqNum = 0;
            UINT32 seqLen = 0;
            UINT64 seqOff = 0;
            for (int i=1; i<line->numFields; ++i)
            {
                if (strncmp(line->fields[i], "SN:", 3) == 0)
                {
                    seqID = line->fields[i]+3;
                    seqNum = count;
                    count += 1;
                }
                else if (strncmp(line->fields[i], "LN:", 3) == 0)
                {
                    seqLen = (UINT32)padLength(str2int(line->fields[i]+3));
                    seqOff = totalLen;
                    totalLen += (UINT64)(seqLen+1);
                }
            }

            // Unless we are marking dups, we don't need to use sequence numbers.
            if (!state->acceptDups)
            {
                // grow seqLens and seqOffs arrays
                if(seqNum % 32768 == 1)
                {
                    state->seqLens = (UINT32*)realloc(state->seqLens, (seqNum+32768)*sizeof(UINT32));
                    state->seqOffs = (UINT64*)realloc(state->seqOffs, (seqNum+32768)*sizeof(UINT64));
                }

                state->seqs[strdup(seqID)] = seqNum;
                state->seqLens[seqNum] = seqLen;
                state->seqOffs[seqNum] = seqOff;
            }
        }
        // Output the header line.
        writeLine(line, state->outputFile);
        if (state->discordantFile != NULL)
            writeLine(line, state->discordantFile);
        if (state->splitterFile != NULL)
            writeLine(line, state->splitterFile);
        disposeSplitLines(line);

        // Read in next line.
        line = readLine(state->inputFile);
    }

    // Output the @PG header lines.
    printPGsamLine(state->outputFile, state);
    printPGsamLine(state->discordantFile, state);
    printPGsamLine(state->splitterFile, state);

    // Make sure we have a header.
    if (count == 1 && !state->acceptDups)
    {
        fatalError("samblaster: Missing header on input sam file.  Exiting.\n");
    }

    // Don't count the "*" entry.
    fprintf(stderr, "samblaster: Loaded %d header sequence entries.\n", count-1);

    // Make sure we have any more lines to process.
    if (line == NULL) return 0;

    // We can now calculate the size of the signatures array, and intialize it.
    if (!state->acceptDups)
    {
        // Make sure we can handle the number of sequences.
        int binCount = (totalLen >> BIN_SHIFT);
        if (binCount >= (1 << 15))
        {
            fatalError("samblaster: Too many sequences in header of input sam file.  Exiting.\n");
        }

        state->binCount = binCount;
        state->sigArraySize = (binCount * 2 + 1) * (binCount * 2 + 1)+1;
#ifdef DEBUG
        fprintf(stderr, "samblaster: Signature Array Size = %d\n", state->sigArraySize);
#endif
        state->sigs = (sigSet_t *) malloc(state->sigArraySize * sizeof(sigSet_t));
        if (state->sigs == NULL) fatalError("samblaster: Unable to allocate signature set array.");
        for (UINT32 i=0; i<state->sigArraySize; i++) hashTableInit(&(state->sigs[i]));
    }

    // Now start processing the alignment records.
    // line already has the first such line in it unless the file only has a header.
    // Keep a ptr to the end of the current list.
    splitLine_t * last = line;
    count = 1;
    while (true)
    {
        splitLine_t * nextLine = readLine(state->inputFile);
        if (nextLine == NULL) break;
        if (strcmp(line->fields[QNAME], nextLine->fields[QNAME]) != 0)
        {
            processSAMBlock(line, state);
            last = line = nextLine;
            count = 1;
        }
        else
        {
            // If we want the SAM lines to be in the output in the same order as the input....
            //     we need to be a little careful about how we make the list.
            last->next = nextLine;
            last = nextLine;
            count += 1;
        }
    }
    // We need to process the final set.
    processSAMBlock(line, state);

    if (!state->quiet)
    {
        if (state->discordantFile != NULL)
            fprintf(stderr, "samblaster: Output %"PRIu64" discordant read pairs to %s\n", discCount/2, state->discordantFileName);
        if (state->splitterFile != NULL)
            fprintf(stderr, "samblaster: Output %"PRIu64" split reads to %s\n", splitCount, state->splitterFileName);
        if (state->unmappedClippedFile != NULL)
            fprintf(stderr, "samblaster: Output %"PRIu64" unmapped/clipped reads to %s\n", unmapClipCount, state->unmappedClippedFileName);
    }

    // Output stats.
    if (state->ignoreUnmated)
    {
        fprintf(stderr, "samblaster: Found %"PRIu64" of %"PRIu64" (%4.2f%%) read ids unmated\n",
                unmatedCount, idCount, ((double)100)*unmatedCount/idCount);
        if (unmatedCount > 0) fprintf(stderr, "samblaster: Please double check that input file is read-id (QNAME) grouped\n");
    }
    if (state->removeDups)
    {
        fprintf(stderr, "samblaster: Removed %"PRIu64" of %"PRIu64" (%4.2f%%) read ids as duplicates",
                dupCount, idCount, ((double)100)*dupCount/idCount);
    }
    else
    {
        fprintf(stderr, "samblaster: Marked %"PRIu64" of %"PRIu64" (%4.2f%%) read ids as duplicates",
                dupCount, idCount, ((double)100)*dupCount/idCount);
    }
    if ((TIMING == 0) || state->quiet)
    {
        fprintf(stderr, ".\n");
    }
    else
    {
#if TIMING
        getrusage(RUSAGE_SELF, &usagebuf);
        time_t endTime = time(NULL);
        struct timeval endRUTime = usagebuf.ru_utime;
        fprintf(stderr, " using %luk memory in ", usagebuf.ru_maxrss);
        fprintTimeMicroSeconds(stderr, diffTVs(&startRUTime, &endRUTime), 3);
        fprintf(stderr, " CPU seconds and ");
        fprintTimeSeconds(stderr, (endTime-startTime), 0);
        fprintf(stderr, " wall time.\n");
#endif
    }

    // Close files.
    if (!streq(state->inputFileName, "stdin")) fclose(state->inputFile);
    if (!streq(state->outputFileName, "stdout")) fclose(state->outputFile);
    if (state->discordantFile != NULL) fclose(state->discordantFile);
    if (state->splitterFile != NULL) fclose(state->splitterFile);
    if (state->unmappedClippedFile != NULL) fclose(state->unmappedClippedFile);

    // Clean up the heap.
    cleanUpSplitLines();
    deleteState(state);
    freeHashTableNodes();
    return 0;
}
