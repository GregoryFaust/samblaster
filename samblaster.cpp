/* -*- Mode: C++ ; indent-tabs-mode: nil ; c-file-style: "stroustrup" -*-

    Project: Fast (but limited) mark duplicates
             Also, optionally pull discordants and splitters.
    Author:  Greg Faust (gf4ea@virginia.edu)
    Date:    October 2013

*/

#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>

// Rename common integer types.
// I like having these shorter name.
typedef uint64_t UINT64;
typedef uint32_t UINT32;

// Some helper routines.
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
#ifndef COMPILE_USER_MODE
#define TIMING
#endif
#ifdef TIMING

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
#endif // TIMING 
#include "Timing.inl"

///////////////////////////////////////////////////////////////////////////////
// Split Lines
///////////////////////////////////////////////////////////////////////////////

// The structure to store "split" input lines, especially SAM lines.
// They form a singly linked list so that we can form groups of them,
//     and also so that we can keep a freelist of them.

// We need to pre-define these for the SAM specific fields.
typedef UINT32 pos_t; // Type for reference offsets.
typedef UINT64 sig_t; // Type for signatures for offsets and lengths.
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
    int   SQO;
    int   EQO;
    bool  discordant;
    bool  splitter;
};

// Creator for splitLine
splitLine_t * makeSplitLine()
{
    splitLine_t * line = (splitLine_t *)malloc(sizeof(splitLine_t));
    line->next = NULL;
    line->bufLen = 0;
    line->maxBufLen = 1000;
    line->buffer = (char *)malloc(line->maxBufLen);
    line->numFields = 0;
    line->maxFields = 100;
    line->fields = (char **)malloc(line->maxFields * sizeof(char *));
    return line;
}

// Destructor for split line.
void disposeSplitLine(splitLine_t * line)
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
        disposeSplitLine(l);
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
    if (splitLineFreeList ==  NULL)
    {
        return makeSplitLine();
    }
    splitLine_t * line = splitLineFreeList;
    splitLineFreeList = splitLineFreeList->next;
    line->next = NULL;
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
        if ((size_t)(line->bufLen + move) > line->maxBufLen)
        {
            fatalError("samblaster: New buffer length exceeds maximum while changing field value.\n");
        }
        // Make the right sized hole for the copy.
        int distance = 1 + line->bufLen - (fp - line->buffer) - move;
        memmove(fp+newLen, fp+oldLen, distance);
        // Correct the total length of the buffer.
        line->bufLen += move;
        // We need to correct the other ptrs as well.
        for (int i=fnum+1; i<line->numFields; i++) line->fields[i] += move;
    }
    // Copy in the new value.
    memcpy(fp, newValue, newLen);
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

///////////////////////////////////////////////////////////////////////////////
// SAM and signature set related structures.
///////////////////////////////////////////////////////////////////////////////

// This controls if we are using hash tables or btrees for the signature structure.
#define useHashTable

#include <map>
#ifdef useHashTable
#include <unordered_set>
#else
#include <set>
#endif

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

// Function needed to get char * map to work.
struct less_str
{
   bool operator()(char const *a, char const *b) const
   {
      return strcmp(a, b) < 0;
   }
};

struct equal_str
{
   bool operator()(char const *a, char const *b) const
   {
      return strcmp(a, b) == 0;
   }
};

// Typedefs for our three main storage structures, the sequence map, the length map, and the signature set.

// We use a map instead of a hash map for sequence names.
// This is because the default hash function on char * hashes the ptr values.
// So, we would need to define our own hash on char * to get things to work properly.
// Not worth it for a structure holding so few members.

// This stores the map between sequence names and sequence numbers.
typedef std::map<const char *, int, less_str> seqMap_t;
// This stores the map between lens and the "length number".
// This is needed as we only have 4 bits to store lengths.
// But most bam files have very few different length reads (unless they do strict quality clipping).
typedef std::map<int, int> lenMap_t;

inline void addSeq(seqMap_t * seqs, char * item, int val)
{
    (*seqs)[item] = val;
}

#ifdef useHashTable
typedef std::unordered_set<sig_t> sigSet_t;
#else
typedef std::set<sig_t> sigSet_t;
#endif
typedef std::pair<sigSet_t::iterator, bool> sigSet_Insert_Return_t;

inline int str2int (char * str)
{
    return strtol(str, NULL, 0);
}

inline pos_t str2pos (char * str)
{
    return strtoul(str, NULL, 0);
}

inline bool checkFlag(splitLine_t * line, int bits)
{
    return (line->flag & bits) != 0;
}

inline void setFlag(splitLine_t * line, int bits)
{
    line->flag |= bits;
}

// Temp buffer to use to form new flag field when marking dups.
char tempBuf[10];
inline void markDup(splitLine_t * line)
{
    setFlag(line, 0x400);
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
    fprintf(output, "_%d\t", checkFlag(line, 0x40) ? 1 : 2);
    outputString(line->fields[1], output);
    fprintf(output, "\n");
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
    sigSet_t *     sigs;
    seqMap_t       seqs;
    splitLine_t ** splitterArray;
    int            splitterArrayMaxSize;
    int            sigArraySize;
    int            minNonOverlap;
    int            maxSplitCount;
    bool           acceptDups;
    bool           excludeDups;
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
    s->sigs = NULL;
    s->minNonOverlap = 20;
    s->maxSplitCount = 2;
    s->acceptDups = false;
    s->excludeDups = false;
    // Used as a temporary location for ptrs to splitter for sort routine.
    s->splitterArrayMaxSize = 1000;
    s->splitterArray = (splitLine_t **)(malloc(s->splitterArrayMaxSize * sizeof(splitLine_t *)));
    return s;
}

void disposeState(state_t * s)
{
    free(s->splitterArray);
    if (s->sigs != NULL) delete[] s->sigs;
    for (seqMap_t::iterator iter = s->seqs.begin(); iter != s->seqs.end(); ++iter) 
    {
        free((char *)(iter->first));
    }
    delete s;
}


///////////////////////////////////////////////////////////////////////////////
// Signatures
///////////////////////////////////////////////////////////////////////////////

inline sig_t calcSig(splitLine_t * first, splitLine_t * second)
{
    // Total nonsense to get the compiler to actually work.
    UINT64 t1 = first->pos;
    UINT64 t2 = t1 << 32;
    UINT64 final = t2 | second->pos;
    return (sig_t)final;
}

inline int calcSigArrOff(splitLine_t * first, splitLine_t * second, seqMap_t & seqs)
{
    int s1 = (first->seqNum * 2) + checkFlag(first, 0x10);
    int s2 = (second->seqNum * 2) + checkFlag(second, 0x10);
    int retval = (s1 * seqs.size() * 2) + s2;
#ifdef DEBUG
    fprintf(stderr, "1st %d %d -> %d 2nd %d %d -> %d count %d result %d\n", 
            firstSeqNum, firstRevBit, s1, secondSeqNum, secondRevBit, s2, numOfSeqs, retval);
#endif
    return retval;
}

///////////////////////////////////////////////////////////////////////////////
// Sequences
///////////////////////////////////////////////////////////////////////////////

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
// Process SAM Blocks
///////////////////////////////////////////////////////////////////////////////

void outputSAMBlock(splitLine_t * block, FILE * output)
{
    for (splitLine_t * line = block; line != NULL; line = line->next)
    {
        writeLine(line, output);
    }
    disposeSplitLines(block);
}

// Some extra fields for in-house statistics.
#ifndef COMPILE_USER_MODE
UINT64 idCount = 0;
UINT64 dupCount = 0;
UINT64 discCount = 0;
UINT64 splitCount = 0;
#endif
// This is the main workhorse that determines if lines are dups or not.
void markDupsDiscordants(splitLine_t * block, state_t * state)
{
    splitLine_t * first = NULL;
    splitLine_t * second = NULL;    
    int count = 0;
    for (splitLine_t * line = block; line != NULL; line = line->next)
    {
        count += 1;
        // Reset discordant flag.
        line->discordant = false;
        // Do this conversion once and store the result.
        line->flag = str2int(line->fields[FLAG]);
        // We don't make our duplicate decisions based on secondaries.
        if (checkFlag(line, 0x100)) continue;
        if      (checkFlag(line, 0x40)) first  = line;
        else if (checkFlag(line, 0x80)) second = line;
    }
    // Make sure we HAVE a first and second of pair.
    if (first == NULL || second == NULL)
    {
        char * temp;
        asprintf(&temp, "samblaster: Can't find first and/or second of pair in sam block of length %d for id: %s\n%s\n",
                 count, block->fields[QNAME], "Are you sure the input is sorted by read ids?");
        fatalError(temp);
    }
    // Never mark pairs as dups if both sides are unmapped.
    if (checkFlag(first, 0x4) && checkFlag(second, 0x4))
    {
        return;
    }
    if (!state->acceptDups)
    {
        // Calculate and store the positions.
        first->pos = str2pos(first->fields[POS]);
        second->pos = str2pos(second->fields[POS]);    
        // Find and store the two seqs numbers.
        // first->seqNum = state->seqs.find(first->fields[RNAME])->second;
        first->seqNum = getSeqNum(first, RNAME, state);
        if (streq(first->fields[RNEXT], "=")) second->seqNum = first->seqNum;
        // else                                 second->seqNum = state->seqs.find(first->fields[RNEXT])->second;
        else                                 second->seqNum = getSeqNum(first, RNEXT, state);
        // The fact of which alignment is first or second in the template is not relevant for determining dups.
        // Therefore, we normalize the pairs based on their characteristics.
        // For orphans, the mapped alignment is first to make sure strand is handled properly.
        // Otherwise, sort by pos, and if equal, sort by sequence num.
        if ((checkFlag(first, 0x4) && !checkFlag(second, 0x4)) ||
            (first->pos > second->pos) ||
            ((first->pos == second->pos) && first->seqNum > second->seqNum))
        {
            splitLine_t * temp = first;
            first = second;
            second = temp;
        }
        // Now find the signature of the pair.
        sig_t sig = calcSig(first, second);
        // Calculate the offset into the signatures array.
        int off = calcSigArrOff(first, second, state->seqs);
        // Attempt insert into the sigs structure.
        // The return value will tell us if it was already there.
        sigSet_Insert_Return_t result = state->sigs[off].insert(sig);
        // Check if the insertion actually happened.
        if (!result.second) 
        {
#ifndef COMPILE_USER_MODE
            dupCount += 1;
#ifdef DEBUG
            fprintf(stderr, "Marking as dup %s\n", line->fields[QNAME]);
#endif
#endif
            // We always mark all or none of a block as dup.
            for (splitLine_t * line = block; line != NULL; line = line->next)
            {
                markDup(line);
            }
        }
    }
    // The first and second help us mark the discordants.
    // Both sides mapped, but pair not properly aligned.
    if (!checkFlag(first, 0x2) && !checkFlag(first, 0x4) && !checkFlag(first, 0x8))
    {
        first->discordant = true;
        second->discordant = true;
    }
}

// Calculate the query offsets from the CIGAR and the strand.
void calcQOs(splitLine_t * line)
{
    int SQO = 0;
    int EQO = 0;
    int QLen = 0;
    int num = 0;
    bool first = true;
    for (int i=0; i<(int)strlen(line->fields[CIGAR]); i++)
    {
        char curchar = line->fields[CIGAR][i];
        int digit = curchar - '0';
        if (digit >= 0 && digit <= 9) num = num * 10 + digit;
        else
        {
            if (curchar == 'H' || curchar == 'S') 
            {
                QLen += num;
                if (first) {SQO += num; EQO += num;}
            }
            else if (curchar == 'M' || curchar == 'I')
            {
                QLen += num;
                EQO += num;
            }
            num = 0;
            first = false;
        }
    }
    // EQO now has a length instead of an offset, so subtract one to make it zero based offset.
    EQO -= 1;
    // Strand normalize the SQO and EQO.
    if (checkFlag(line, 0x10))
    {
        int temp = EQO;
        EQO = QLen - (SQO + 1);
        SQO = QLen - (temp + 1);
    }
    line->SQO = SQO;
    line->EQO = EQO;
}

// Sort ascending in SQO.
int compQOs(const void * p1, const void * p2)
{
    splitLine_t * l1 = (*(splitLine_t **)p1);
    splitLine_t * l2 = (*(splitLine_t **)p2);
    return (l1->SQO - l2->SQO);
}

void markSplitter(splitLine_t * block, state_t * state, int read)
{
    // Count the secondaries we have for this read (if any), and store their ptrs into an array.
    int count = 0;
    for (splitLine_t * line = block; line != NULL; line = line->next)
    {
        // For all the ones that are the current read of interest....
        if (checkFlag(line, read))
        {
            // Mark the line as not a splitter.
            line->splitter = false;
            // If it is not a secondary.
            if (!checkFlag(line, 0x100))
            {
                // Make sure the primary is mapped!
                if (checkFlag(line, 0x4)) return;
            }
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
    if (count < 2 || count > state->maxSplitCount) return;
    // We have the lines of interest in an array.
    // Calculate the query positions.
    for (int i=0; i<count; i++)
    {
        calcQOs(state->splitterArray[i]);
    }
    // We need to sort it by strand normalized query offsets.
    qsort(state->splitterArray, count, sizeof(splitLine_t *), compQOs);
    // Now check for the amount of overlap.
    splitLine_t * l1 = state->splitterArray[0];
    splitLine_t * l2;
    for (int i=1; i<count; i++)
    {
        l2 = state->splitterArray[i];
        // We don't allow negative overlap, as that will lead to wrong non-overlap calculation.
        int overlap = std::max(1 + std::min(l1->EQO, l2->EQO) - std::max(l1->SQO, l2->SQO), 0);
        int alen1 = 1 + l1->EQO - l1->SQO;
        int alen2 = 1 + l2->EQO - l2->SQO;
        int mno = std::min(alen1-overlap, alen2-overlap);
        if (mno >= state->minNonOverlap) 
        {
            l1->splitter = true;
            l2->splitter = true;
        }

#ifdef DEBUG
        fprintf(stderr, "After sort, for l1, QOs are %d and %d\n", l1->SQO, l1->EQO);
        fprintf(stderr, "After sort, for l2, QOs are %d and %d\n", l2->SQO, l2->EQO);
        fprintf(stderr, "After sort, overlap = %d,  mno = %d l1split %d l2split %d\n", overlap, mno, l1->splitter, l2->splitter);
#endif

        // Get ready for the next iteration.
        l1 = l2;
    }
}

void processSAMBlock(splitLine_t * block, state_t * state)
{
#ifndef COMPILE_USER_MODE
    idCount += 1;
#endif
    // First mark dups and find the discordants.
    // These share a lot of looking at flag bits, so make sense to put together.
    markDupsDiscordants(block, state);
    // Look for splitters.
    // Since this is expensive, don't do it unless the user asked for them.
    if (state->splitterFile != NULL)
    {
        // Check the first read for splitter.
        markSplitter(block, state, 0x40);
        // Check the second read for splitter.
        markSplitter(block, state, 0x80);
    }
    // Now do the output.
    for (splitLine_t * line = block; line != NULL; line = line->next)
    {
        // Write to the output file.
        writeLine(line, state->outputFile);
        // Write to discordant file if appropriate.
        if (state->discordantFile != NULL && line->discordant && !(state->excludeDups && checkFlag(line, 0x400)))
        {
            writeLine(line, state->discordantFile);
#ifndef COMPILE_USER_MODE
            discCount += 1;
#endif
        }
        // Write to splitter file if appropriate.
        if (state->splitterFile != NULL && line->splitter && !(state->excludeDups && checkFlag(line, 0x400)))
        {
            writeSAMlineWithIdNum(line, state->splitterFile);
#ifndef COMPILE_USER_MODE
            splitCount += 1;
#endif
        }
    }
    disposeSplitLines(block);
}


///////////////////////////////////////////////////////////////////////////////
// Main Routine with helpers.
///////////////////////////////////////////////////////////////////////////////

void printVersionString()
{
    fprintf(stderr, "samblaster: Version 0.1.%d\n", BUILDNUM);
}

void printUsageString()
{
    const char* useString = "Tool to mark duplicates and optionally output split reads and/or discordant pairs.\n"
        "Input sam file must contain paired end data, contain sequence header and be sorted by read ids.\n"
        "Output will be all alignments in the same order as input, with duplicates marked with FLAG 0x400.\n\n"

        "Usage:\n"
        "For use with as a post process on an aligner (eg. bwa mem):\n"
        "     bwa mem index samp.r1.fq sam.r2.fq | samblaster [-s samp.split.sam] [-d samp.disc.sam] | samtools view -Sb - > samp.out.bam\n"
        "For use with a pre-existing bam file to pull split reads and/or discordants:\n"
        "     samtools view -h samp.bam | samblaster [-s samp.split.sam] [-d samp.disc.sam] -o /dev/null\n\n"

        "Input/Output Options:\n"
        "-i --input          FILE Input sam file [stdin].\n"
        "-o --output         FILE Output sam file for all input alignments [stdout].\n"
        "-d --discordantFile FILE Output discordant read pairs to this file. By default, no discordant file is produced.\n"
        "-s --splitterFile   FILE Output split reads to this file abiding by paramters below. By default, no splitter file is produced.\n\n"

        "Other Options:\n"
        "-a --acceptDupMarks      Accept duplicate marks already in input file instead of looking for duplicates in the input.\n"
        "-e --excludeDups         Exclude reads marked as duplicates from discordant and/or splitter file.\n"
        "-c --maxSplitCount  INT  Maximum number of split alignments for a read to be included in splitter file.[2]\n"
        "-m --minNonOverlap  INT  Minimum non-overlaping base pairs between two alignments for a read to be included in splitter file. [20]\n";

        printVersionString();
        fprintf(stderr, useString);
}

void printUsageStringAbort()
{
    printUsageString();
    exit(1);
}

int main (int argc, char *argv[])
{
    // Set up timing.
#ifdef TIMING
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
        else if (streq(argv[argi], "-c") || streq(argv[argi],"--maxSplitCount"))
        {
            argi++;
            state->maxSplitCount = str2int(argv[argi]);
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
        else if (streq(argv[argi], "-m") || streq(argv[argi],"--minNonOverlap"))
        {
            argi++;
            state->minNonOverlap = str2int(argv[argi]);
        }
        else if (streq(argv[argi], "-o") || streq(argv[argi],"--output"))
        {
            argi++;
            state->outputFileName = argv[argi];
        }
        else if (streq(argv[argi], "-s") || streq(argv[argi],"--splitterFile"))
        {
            argi++;
            state->splitterFileName = argv[argi];
        }
        else
        {
            fprintf(stderr, "samblaster: Unrecognized option: %s\n", argv[argi]);
            printUsageString();
            return 1;
        }
    }

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
    if (!streq(state->inputFileName, "stdin"))
    {
        fprintf(stderr, "samblaster: Opening %s for read.\n", state->inputFileName);
        state->inputFile = fopen(state->inputFileName, "r");
        if (state->inputFile == NULL) fsError(state->inputFileName);
    }
    if (!streq(state->outputFileName, "stdout"))
    {
        fprintf(stderr, "samblaster: Opening %s for write.\n", state->outputFileName);
        state->outputFile = fopen(state->outputFileName, "w");
        if (state->outputFile == NULL) fsError(state->outputFileName);
    }
    if (!streq(state->discordantFileName, ""))
    {
        fprintf(stderr, "samblaster: Opening %s for write.\n", state->discordantFileName);
        state->discordantFile = fopen(state->discordantFileName, "w");
        if (state->discordantFile == NULL) fsError(state->discordantFileName);
    }
    if (!streq(state->splitterFileName, ""))
    {
        fprintf(stderr, "samblaster: Opening %s for write.\n", state->splitterFileName);
        state->splitterFile = fopen(state->splitterFileName, "w");
        if (state->splitterFile == NULL) fsError(state->splitterFileName);
    }
    fprintf(stderr, "samblaster: Inputting from %s and outputting to %s\n", state->inputFileName, state->outputFileName);

    // Read in the SAM header and create the seqs structure.
    state->seqs[strdup("*")] = 0;
    int count = 1;
    splitLine_t * line;
    while (true)
    {
        line = readLine(state->inputFile);
        // Check if we have exhausted the header.
        if (line == NULL || line->fields[0][0] != '@') break;
        // Process input line to see if it defines a sequence.
        if (streq(line->fields[0], "@SQ"))
        {
            for (int i=1; i<line->numFields; ++i)
            {
                if (strncmp(line->fields[i], "SN:", 3) == 0)
                {
                    state->seqs[strdup(line->fields[i]+3)] = count;
                    count += 1;
                }
            }
        }
        // Output the header line.
        writeLine(line, state->outputFile);
        if (state->discordantFile != NULL)
            writeLine(line, state->discordantFile);
        if (state->splitterFile != NULL)
            writeLine(line, state->splitterFile);
        disposeSplitLines(line);
    }
    
    // Make sure we have a header.
    if (count == 1 && !state->acceptDups)
    {
        fatalError("samblaster: Missing header on input sam file.  Exiting.\n");
    }

#ifndef COMPILE_USER_MODE
    // Don't count the "*" entry.
    fprintf(stderr, "samblaster: Loaded %d header sequence entries.\n", count-1);
#endif

    // Make sure we have any more lines to process.
    if (line == NULL) return 0;

    // We can now calculate the size of the signatures array, and intialize it.
    if (!state->acceptDups)
    {
        state->sigArraySize = count * count * 4;
        state->sigs = new sigSet_t[state->sigArraySize];
    }

    // Now start processing the alignment records.
    // line already has the first such line in it unless the file only has a header.
    // Initialize our signature set.
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

#ifndef COMPILE_USER_MODE
    if (state->discordantFile != NULL)
        fprintf(stderr, "samblaster: Output %lu discordant pairs to %s\n", discCount/2, state->discordantFileName);
    if (state->splitterFile != NULL)
        fprintf(stderr, "samblaster: Output %lu split reads to %s\n", splitCount, state->splitterFileName);
#ifdef TIMING
    // Output stats.
    time_t endTime = time(NULL);
    getrusage(RUSAGE_SELF, &usagebuf);
    struct timeval endRUTime = usagebuf.ru_utime;    
    fprintf(stderr, "samblaster: Marked %lu of %lu (%4.2f%%) read ids as duplicates", dupCount, idCount, ((double)100)*dupCount/idCount);
    fprintf(stderr, " using %luk memory in ", usagebuf.ru_maxrss);
    fprintTimeMicroSeconds(stderr, diffTVs(&startRUTime, &endRUTime), 3);
    fprintf(stderr, " CPU seconds and ");
    fprintTimeSeconds(stderr, (endTime-startTime), 0);
    fprintf(stderr, " wall time.\n");
#else
    struct rusage usagebuf;
    getrusage(RUSAGE_SELF, &usagebuf);
    fprintf(stderr, "samblaster: Marked %lu of %lu (%4.2f%%) read ids as duplicates", dupCount, idCount, ((double)100)*dupCount/idCount);
    fprintf(stderr, " using %luk memory.\n", usagebuf.ru_maxrss);
#endif // TIMING
#endif // COMPILE_USER_MODE

    // Close files.
    if (!streq(state->inputFileName, "stdin")) fclose(state->inputFile);
    if (!streq(state->outputFileName, "stdout")) fclose(state->outputFile);
    if (state->discordantFile != NULL) fclose(state->discordantFile);
    if (state->splitterFile != NULL) fclose(state->splitterFile);
    // These clean everything up, but they take forever to run.
    // So, just punt and leave the mess behind.
#ifdef DEBUG
    cleanUpSplitLines();
    disposeState(state);
#endif
    return 0;
}
