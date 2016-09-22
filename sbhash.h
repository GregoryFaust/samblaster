/* -*- Mode: C++ ; indent-tabs-mode: nil ; c-file-style: "stroustrup" -*-

    Project: samblaster
             Fast mark duplicates in read-ID grouped SAM file.
             Also, optionally pull discordants, splitters, and/or unmappend/clipped reads.
    Author:  Greg Faust (gf4ea@virginia.edu)
    Date:    October 2013

    File:    sbhash.h  header file for our hash table.

    License Information:

    Copyright 2013-2016 Gregory G. Faust

    Licensed under the MIT license (the "License");
    You may not use this file except in compliance with the License.
    You may obtain a copy of the License at http://opensource.org/licenses/MIT

*/

#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>

// Rename common integer types.
// I like having these shorter name.
typedef uint64_t UINT64;
typedef uint32_t UINT32;

///////////////////////////////////////////////////////////////////////////////
// Hash Table Collision Nodes
///////////////////////////////////////////////////////////////////////////////

// This controls the number of data values to store in each hadnNode.
#define HASHNODE_PAYLOAD_SIZE 3

typedef struct hashNode hashNode_t;
struct hashNode
{
    hashNode_t * next;
    UINT64 values[HASHNODE_PAYLOAD_SIZE];
};

hashNode_t * getHashNode();
void disposeHashNode(hashNode_t * node);

///////////////////////////////////////////////////////////////////////////////
// Hash Table
///////////////////////////////////////////////////////////////////////////////

typedef struct hashTable hashTable_t;
struct hashTable
{
    UINT64 * table;
    UINT32   size;
    UINT32   entries;
    ~hashTable();
};

hashTable_t * makeHashTable();
void deleteHashTable(hashTable_t * ht);
bool hashTableInsert(hashTable_t * ht, UINT64 value);
void hashTableInit(hashTable_t * ht, int size=0);
void freeHashTableNodes();
