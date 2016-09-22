/* -*- Mode: C++ ; indent-tabs-mode: nil ; c-file-style: "stroustrup" -*-

    Project: samblaster
             Fast mark duplicates in read-ID grouped SAM file.
             Also, optionally pull discordants, splitters, and/or unmappend/clipped reads.
    Author:  Greg Faust (gf4ea@virginia.edu)
    Date:    October 2013

    File:    sbhash.cpp  code file for our hash table.

    License Information:

    Copyright 2013-2016 Gregory G. Faust

    Licensed under the MIT license (the "License");
    You may not use this file except in compliance with the License.
    You may obtain a copy of the License at http://opensource.org/licenses/MIT

*/

#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <sys/mman.h>
#include <cmath>
#include "sbhash.h"

///////////////////////////////////////////////////////////////////////////////
//     Slab Allocator common code
///////////////////////////////////////////////////////////////////////////////
//     We will lazily allocate more slabs when needed, and only clean up at the end.
//////////////////////////////////////////////////////////////////////////////

void fatalError(const char * errorStr);
void checkFSerrWithFilename (ssize_t returnCode)
{
    if (returnCode == -1)
    {
        char * temp;
        if (errno == ENOMEM)
            temp = (char *)"samblaster: Insufficient memory available to satisfy allocation request.\n";
        else
            asprintf(&temp, "File system error %d trying to allocate or free memory\n", errno);
        fatalError(temp);
    }
}

// Allocate big blocks of memory.
char * blockMalloc(ssize_t size)
{
    char * retval = (char *)mmap(NULL, size, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANON, 0, 0);
    checkFSerrWithFilename((ssize_t)retval);
    return retval;
}

// Free big blocks of memory.  Not the size is needed.
void blockFree(char * ptr, ssize_t size)
{
    int err = munmap(ptr, size);
    checkFSerrWithFilename(err);
}

typedef struct LBMallocBlock LBMallocBlock_t;
struct LBMallocBlock
{
    char *            block;     // Pointer to the payload block
    LBMallocBlock_t * next;      // Pointer to the next allocation block
    size_t            size;      // Size of the allocated block
};

char * pushNewLBMallocBlock(int blockSize, LBMallocBlock_t **blockArrayPtr)
{
    char * newBlock = blockMalloc(blockSize);
    LBMallocBlock_t * newMallocBlock = (LBMallocBlock_t *)malloc(sizeof(LBMallocBlock_t));
    if (newMallocBlock == NULL) fatalError("samblaster: Insufficeint memory available to allocate (more) objects.");
    newMallocBlock->size = blockSize;
    newMallocBlock->block = newBlock;
    newMallocBlock->next = *blockArrayPtr;
    *blockArrayPtr = newMallocBlock;
    return newBlock;
}

void freeLBMallocBlocks(LBMallocBlock_t * block)
{
    while (block != NULL)
    {
        LBMallocBlock_t * nextBlock = block->next;
        blockFree(block->block, block->size);
        free(block);
        block = nextBlock;
    }
}

///////////////////////////////////////////////////////////////////////////////
// Hash Table Collision Nodes
///////////////////////////////////////////////////////////////////////////////

#define newNodeCount 4096

// Ptr to head of linked list of allocated node slabs.
LBMallocBlock_t * nodeBlockList = NULL;
// Ptr to head of linked list of free node objects.
hashNode_t * hashNodeFreeList = NULL;

void makeMoreHashNodes()
{
    hashNode_t * nodeArray = (hashNode_t *)pushNewLBMallocBlock(sizeof(hashNode_t) * newNodeCount, &nodeBlockList);
    for (int i=1; i<newNodeCount; i++)
    {
        (nodeArray + (i - 1))->next = (nodeArray + i);
    }
    (nodeArray + (newNodeCount - 1))->next = NULL;
    hashNodeFreeList = nodeArray;
}

hashNode_t * getHashNode()
{
    if (hashNodeFreeList == NULL) makeMoreHashNodes();
    hashNode_t * node = hashNodeFreeList;
    hashNodeFreeList = hashNodeFreeList->next;
    node->next = NULL;
    for (int i=0; i<HASHNODE_PAYLOAD_SIZE; i++) node->values[i] = 0;
    return node;
}

// I don't think this is currently being called, as we always put the entire string of nodes on the freelist.
void disposeHashNode(hashNode_t * node)
{
    node->next = hashNodeFreeList;
    hashNodeFreeList = node;
}

void freeHashTableNodes()
{
    freeLBMallocBlocks(nodeBlockList);
}

///////////////////////////////////////////////////////////////////////////////
// Hash Table
///////////////////////////////////////////////////////////////////////////////
// We are going to depend on an old hack.
// ptrs to 8 byte things will be 8 byte aligned.
// Therefore, the lower 3 bits will be zero.
// Also, no known chromosome offset requires all 32 bits.
// So we will roll the signature up one bit and put a one there.
// We can then tell apart the three following state for a table entry:
//  0         -> empty bucket
//  low bit 1 -> value
//  low bit 0 -> ptr to overflow nodes.
///////////////////////////////////////////////////////////////////////////////

inline UINT64 makeValue(UINT64 value)
{
    return (value << 1) | 1;
}

inline UINT64 unmakeValue(UINT64 value)
{
    return (value >> 1);
}

inline hashNode_t * makePtr(UINT64 value)
{
    return (hashNode_t *)value;
}

inline bool isEmpty(UINT64 value)
{
    return (value == 0);
}

inline bool isValue(UINT64 value)
{
    return ((value & 1) != 0);
}

#define numOfSizes 27
static UINT32 hashTableSizes [] = {0, 23, 47, 97, 199, 409, 823, 1741, 3739, 7517, 15173, 30727, 62233, 126271, 256279, 520241, 1056323,
                                   2144977, 4355707, 8844859, 17961079, 36473443, 74066549, 150406843, 305431229, 620239453, 1259520799};

inline UINT32 hash(UINT64 value)
{
    return (UINT32)value;
}

void hashTableInit(hashTable_t * ht, int size)
{
    ht->entries = 0;
    ht->size = size;
    if (size == 0)
    {
        ht->table = (UINT64 *)NULL;
        return;
    }
    ht->table = (UINT64 *)calloc(ht->size, sizeof(UINT64));
    if (ht->table == NULL) fatalError("samblaster: unable to allocate hash table.\n");
}

hashTable_t * makeHashTable()
{
    hashTable_t * ht = (hashTable_t *)malloc(sizeof(hashTable_t));
    if (ht == NULL) fatalError("samblaster: unable to allocate hash table.\n");
    hashTableInit(ht, hashTableSizes[0]);
    return ht;
}

// Use a C++ style destructor so that arrays of hash tables will be cleaned up automagically.
hashTable::~hashTable()
{
    if (table != NULL) free(table);
}

// C style delete.
void deleteHashTable(hashTable_t * ht)
{
    if (ht->table != NULL) free(ht->table);
}

void resizeHashTable(hashTable_t * ht)
{
    // Find out what size table is next.
    int newsize = 0;
    for (int i=0; i<numOfSizes; i++)
    {
        if (hashTableSizes[i] == ht->size)
        {
            newsize = hashTableSizes[i+1];
            break;
        }
    }

    // Remember the current values array.
    UINT64 * oldtable = ht->table;
    int size = ht->size;
    // Now reinit the hash table with a new table, etc.
    hashTableInit(ht, newsize);
    // Now iterate over all values and rehash them into the new table.
    for (int i=0; i<size; i++)
    {
        UINT64 value = oldtable[i];
        if (isEmpty(value)) continue;
        if (isValue(value)) {hashTableInsert(ht, unmakeValue(value)); continue;}
        // We need to iterate through the nodes.
        hashNode_t * node = makePtr(value);
        while (true)
        {
            for (int j=0; j<HASHNODE_PAYLOAD_SIZE; j++)
            {
                value = node->values[j];
                if (isEmpty(value)) break;
                hashTableInsert(ht, unmakeValue(value));
            }
            if (node->next == NULL) break;
            node = node->next;
        }
        // We need to free up the nodes.
        // TODO move out of line.
        node->next = hashNodeFreeList;
        hashNodeFreeList = makePtr(oldtable[i]);
    }

    // Free up the oldtable.
    if (oldtable != NULL) free(oldtable);
}

bool hashTableInsert(hashTable_t * ht, UINT64 value)
{
    // See if we have reached our size limit.
    if (ht->entries == ht->size) resizeHashTable(ht);
    int bucket = hash(value) % ht->size;
    // We need to empty the low order bit so that we can tell the difference between values and ptrs.
    value = makeValue(value);
    UINT64 curvalue = ht->table[bucket];
    // The empty case should be most common.
    if (isEmpty(curvalue))
    {
        ht->table[bucket] = value;
        ht->entries += 1;
        return true;
    }
    // The value case should be next most common.
    if (isValue(curvalue))
    {
        // The value is already here.
        if (curvalue == value) return false;
        // We have a collision and need to add an overflow node.
        hashNode_t * node = getHashNode();
        ht->table[bucket] = (UINT64)node;
        node->values[0] = curvalue;
        // Note that this test doesn't cost us anything as it happens at compile time.
        if (HASHNODE_PAYLOAD_SIZE >= 2)
        {
            node->values[1] = value;
        }
        else
        {
            // We need to add a second new node.
            hashNode_t * secondNode = getHashNode();
            node->next = secondNode;
            secondNode->values[0] = value;
        }
        ht->entries += 1;
        return true;
    }
    // The overflow node case.
    hashNode_t * curNode = makePtr(curvalue);
    while (true)
    {
        for (int i=0; i<HASHNODE_PAYLOAD_SIZE; i++)
        {
            // Check if we have an empty slot.
            if (curNode->values[i] == 0)
            {
                curNode->values[i] = value;
                ht->entries += 1;
                return true;
            }
            // Check if the value matches the current value.
            if (curNode->values[i] == value) return false;
        }
        if (curNode->next == NULL) break;
        curNode = curNode->next;
    }
    // If we are here, we need a new node.
    hashNode_t * node = getHashNode();
    curNode->next = node;
    node->values[0] = value;
    ht->entries += 1;
    return true;
}
