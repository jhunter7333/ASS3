#include "radix.h"
#include "utils.h"
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>

#define INITIALCAPACITY 1024

/* Packed radix tree. */
struct radixTree {
    // Info to derive bits per element.
    int numPieces;
    int height;
    int width;

    // Structure of tree packed as struct of (bit) arrays.
    int nodeCount;
    int nodeCapacity;
 
    int *prefixBitStartBytes;
    
    int *prefixBitsBytes;
    // Bits are stored in one contiguous bit array referred to by all nodes.
    int prefixBitsAllocated;
    int prefixBitsUsed;
    unsigned char *prefixBytes;
    int *branchABytes;
    int *branchBBytes;
};

// Helper structure.
struct radixTreeNode;

struct radixTreeNode {
    int nodeIdx;
    int bitStart;
    int numBits;
    int branchA;
    int branchB;
};

/* Number of bits in a single character. */
#define BITS_PER_BYTE 8
#define SIZE 8U

int getBit(unsigned char *s, unsigned int bitIndex){
    /* bitIndex >= 0 is forced by type. */
    // assert(s && bitIndex >= 0);
    assert(s);
    unsigned int byte = bitIndex / BITS_PER_BYTE;
    unsigned int indexFromLeft = bitIndex % BITS_PER_BYTE;
    /* 
        Since we split from the highest order bit first, the bit we are interested
        will be the highest order bit, rather than a bit that occurs at the end of the
        number. 
    */
    unsigned int offset = (BITS_PER_BYTE - (indexFromLeft) - 1) % BITS_PER_BYTE;
    unsigned char byteOfInterest = s[byte];
    unsigned int offsetMask = (1 << offset);
    unsigned int maskedByte = (byteOfInterest & offsetMask);
    /*
        The masked byte will still have the bit in its original position, to return
        either 0 or 1, we need to move the bit to the lowest order bit in the number.
    */
    unsigned int bitOnly = maskedByte >> offset;
    return bitOnly;
}

int queryRadixMemoryUsage(struct radixTree *tree) {
    int memoryUsage = 0;
    /* 
        Data for each node:
        prefixBitStartByte, prefixBits, branchA, branchB
     */
    memoryUsage += tree->nodeCount * (sizeof(int) * (1 + 1 + 1 + 1));
    /*
        Data used in bits - a bit spilling over one byte takes one more byte.
    */
    memoryUsage += (tree->prefixBitsUsed + (BITS_PER_BYTE - 1)) / BITS_PER_BYTE;

    /*
        Radix tree metadata is not counted - whether this is counted or not is
        fairly reasonable, but is chosen not to be counted here.
        All other data allocated in structure is optional.
    */
    return memoryUsage;
}

struct radixTree *getNewRadixTree(int numPieces, int height, int width) {
    struct radixTree *rt = (struct radixTree *) malloc(sizeof(struct radixTree));
    assert(rt);

    rt->numPieces = numPieces;
    rt->height = height;
    rt->width = width;
    
    rt->nodeCount = 0;
    rt->nodeCapacity = 0;

    rt->prefixBitStartBytes = NULL;
    rt->prefixBitsBytes = NULL;
    rt->prefixBytes = NULL;
    rt->branchABytes = NULL;
    rt->branchBBytes = NULL;

    rt->prefixBitsAllocated = 0;
    rt->prefixBitsUsed = 0;

    return rt;
}

int calcBits(int x) {
    int bitsNeeded;
    /* Highest stored value. */
    x = x - 1;
    for(bitsNeeded = 0; x > 0; bitsNeeded++) {
        x >>= 1;
    }
    return bitsNeeded;
}

struct radixTreeNode getTreeNode(struct radixTree *tree, int idx);

struct radixTreeNode getTreeNode(struct radixTree *tree, int idx){
    struct radixTreeNode rt;
    
    rt.nodeIdx = idx;
    rt.bitStart = tree->prefixBitStartBytes[idx];
    rt.numBits = tree->prefixBitsBytes[idx];
    rt.branchA = tree->branchABytes[idx];
    rt.branchB = tree->branchBBytes[idx];

    return rt;
}

void storeNode(struct radixTree *tree, struct radixTreeNode *node);

void storeNode(struct radixTree *tree, struct radixTreeNode *node){
    if(node->nodeIdx >= tree->nodeCount) {
        while(node->nodeIdx + 1 > tree->nodeCapacity){
            /* Resize needed. */
            tree->nodeCapacity *= 2;
            tree->prefixBitStartBytes = (int *) realloc(tree->prefixBitStartBytes, tree->nodeCapacity * sizeof(int));
            assert(tree->prefixBitStartBytes);
            tree->prefixBitsBytes = (int *) realloc(tree->prefixBitsBytes, tree->nodeCapacity * sizeof(int));
            assert(tree->prefixBitsBytes);
            tree->branchABytes = (int *) realloc(tree->branchABytes, tree->nodeCapacity * sizeof(int));
            assert(tree->branchABytes);
            tree->branchBBytes = (int *) realloc(tree->branchBBytes, tree->nodeCapacity * sizeof(int));
            assert(tree->branchBBytes);
        }
        tree->nodeCount = node->nodeIdx + 1;
    }
    tree->prefixBitStartBytes[node->nodeIdx] = node->bitStart;
    tree->prefixBitsBytes[node->nodeIdx] = node->numBits;
    tree->branchABytes[node->nodeIdx] = node->branchA;
    tree->branchBBytes[node->nodeIdx] = node->branchB;
}

/* Checks if the state is present in the radix tree. */
int checkPresent(struct radixTree *tree, unsigned char *bitPacked, int atomCount) {
    int pBits = calcBits(tree->numPieces);
    int hBits = calcBits(tree->height);
    int wBits = calcBits(tree->width);
    int atomSize = pBits + hBits + wBits;

    /* Full check, so bits contain location of all pieces. */
    int bitCount = atomSize * atomCount;

    if(tree->nodeCount == 0){
        return NOTPRESENT;
    }

    /* Get root. */
    struct radixTreeNode node = getTreeNode(tree, 0);

    /* Search */
    int progress = 0;
    for(int i = 0; i < bitCount; i++) {
        if(progress == node.numBits){
            /* Branch. */
            if(getBit(bitPacked, i) == 0) {
                node = getTreeNode(tree, node.branchA);
                progress = 0;
            } else {
                node = getTreeNode(tree, node.branchB);
                progress = 0;
            }
        }
        if(getBit(tree->prefixBytes, node.bitStart + progress) != getBit(bitPacked, i)){
            /* Mismatch, not in tree. */
            return NOTPRESENT;
        }
        progress++;
    }
    /* Got through whole bitPacked representation. Should be true since we assume bitPacked items are always inserted. */
    assert(progress == node.numBits);
    return PRESENT;
}

/* Write bitCount bits, starting from startBit from the bitPacked value into 
    the prefixBytes of the tree */
void writeNewBits(struct radixTree *tree, unsigned char *bitPacked, int startBit, int bitCount);

// BitOn and BitOff written by Danielle Jayanthy for COMP20007 Semester 1 2025 Assignment 2
// Modified to work consistently with the bit ordering in the getBit function.

// Set the i-th bit to true
void bitOn( unsigned char A[], unsigned int bitIndex ) {
    A[bitIndex/SIZE] |= 1 << (SIZE - (bitIndex%SIZE) - 1);
}

// set the i-th bit to false
void bitOff( unsigned char A[], unsigned int bitIndex ) {
    A[bitIndex/SIZE] &= ~(1 << (SIZE - (bitIndex%SIZE) - 1));
}

void writeNewBits(struct radixTree *tree, unsigned char *bitPacked, int startBit, int bitCount){
    /* Check space for bits. */
    while(tree->prefixBitsUsed + bitCount > tree->prefixBitsAllocated) {
        tree->prefixBitsAllocated *= 2;
        tree->prefixBytes = (unsigned char *) realloc(tree->prefixBytes, tree->prefixBitsAllocated/BITS_PER_BYTE * sizeof(unsigned char));
        assert(tree->prefixBytes);
    }
    for(int i = 0; i < bitCount; i++) {
        if(getBit(bitPacked, startBit + i) == 0){
            bitOff(tree->prefixBytes, tree->prefixBitsUsed + i);
        } else {
            bitOn(tree->prefixBytes, tree->prefixBitsUsed + i);
        }
    }
    tree->prefixBitsUsed += bitCount;
}

/* Inserts the state into the radix tree. */
void insertRadixTree(struct radixTree *tree, unsigned char *bitPacked, int atomCount) {
    int pBits = calcBits(tree->numPieces);
    int hBits = calcBits(tree->height);
    int wBits = calcBits(tree->width);
    int atomSize = pBits + hBits + wBits;

    /* Full check, so bits contain location of all pieces. */
    int bitCount = atomSize * atomCount;

    /* Simple check, do not insert if already present. */
    if(checkPresent(tree, bitPacked, atomCount)){
        return;
    }

    /* Empty tree. */
    if(tree->nodeCapacity == 0) {
        tree->nodeCapacity = INITIALCAPACITY;
        
        tree->prefixBitStartBytes = (int *) malloc(sizeof(int) * INITIALCAPACITY);
        assert(tree->prefixBitStartBytes);
        tree->prefixBitsBytes = (int *) malloc(sizeof(int) * INITIALCAPACITY);
        assert(tree->prefixBitsBytes);
        tree->branchABytes = (int *) malloc(sizeof(int) * INITIALCAPACITY);
        assert(tree->branchABytes);
        tree->branchBBytes = (int *) malloc(sizeof(int) * INITIALCAPACITY);
        assert(tree->branchBBytes);

        (tree->prefixBitStartBytes)[0] = 0;
        (tree->branchABytes)[0] = NOCHILD;
        (tree->branchBBytes)[0] = NOCHILD;
        (tree->prefixBitsBytes)[0] = bitCount;
        
        // Start with INITIALCAPACITY full length prefixes of bitCount length.
        tree->prefixBitsAllocated = ((((bitCount * INITIALCAPACITY + (BITS_PER_BYTE - 1)))/BITS_PER_BYTE)) * sizeof(unsigned char) * BITS_PER_BYTE;
        tree->prefixBytes = (unsigned char *) calloc((((bitCount * INITIALCAPACITY + (BITS_PER_BYTE - 1)))/BITS_PER_BYTE), sizeof(unsigned char));
        assert(tree->prefixBytes);

        writeNewBits(tree, bitPacked, 0, bitCount);

        (tree->nodeCount)++;
        return;
    }

    /* Find mismatch */
    /* Get root. */
    struct radixTreeNode node = getTreeNode(tree, 0);

    /* Search */
    int progress = 0;
    for(int i = 0; i < bitCount; i++) {
        if(progress == node.numBits){
            /* Branch. */
            if(getBit(bitPacked, i) == 0) {
                node = getTreeNode(tree, node.branchA);
                progress = 0;
            } else {
                node = getTreeNode(tree, node.branchB);
                progress = 0;
            }
        }
        if(getBit(tree->prefixBytes, node.bitStart + progress) != getBit(bitPacked, i)){
            /* Mismatch, not in tree. Add to tree. */
            /* Part 0: Root node changes. */
            struct radixTreeNode newRoot;
            newRoot = node;
            newRoot.numBits = progress;
            // newRoot.branchA 
            // newRoot.branchB
            /* Part 1: Node generated from bit packed insertion. */
            struct radixTreeNode newNode;
            newNode.nodeIdx = tree->nodeCount + 1;
            newNode.bitStart = tree->prefixBitsUsed;
            int remainingBits = bitCount - i;
            newNode.numBits = remainingBits;
            newNode.branchA = NOCHILD;
            newNode.branchB = NOCHILD;
            /* Part 2: Node generated from existing tree. */
            struct radixTreeNode existingNode;
            existingNode.nodeIdx = tree->nodeCount;
            existingNode.bitStart = node.bitStart + progress;
            int existingRemaining = node.numBits - progress;
            existingNode.numBits = existingRemaining;
            /* Inherits existing children. */
            existingNode.branchA = node.branchA;
            existingNode.branchB = node.branchB;
            if(getBit(bitPacked, i) == 0) {
                newRoot.branchA = newNode.nodeIdx;
                newRoot.branchB = existingNode.nodeIdx;
            } else {
                newRoot.branchA = existingNode.nodeIdx;
                newRoot.branchB = newNode.nodeIdx;
            }
            /* Store prefix data. */
            writeNewBits(tree, bitPacked, i, remainingBits);
            /* Store nodes. */
            storeNode(tree, &existingNode);
            storeNode(tree, &newNode);
            storeNode(tree, &newRoot);
            return;
        }
        progress++;
    }
}

void writeNewBitsnCr(unsigned char *destBits, int destFilledBits, unsigned char *bitPacked, int startBit, int bitCount);

void writeNewBitsnCr(unsigned char *destBits, int destFilledBits, unsigned char *bitPacked, int startBit, int bitCount){
    for(int i = 0; i < bitCount; i++) {
        if(getBit(bitPacked, startBit + i) == 0){
            bitOff(destBits, destFilledBits + i);
        } else {
            bitOn(destBits, destFilledBits + i);
        }
    }
}

/* =========================
   Lifted helpers (no nesting)
   ========================= */

/* Former inner: packPartial from checkPresentnCr */
static int packPartialCheck(struct radixTree *tree,
                            unsigned char *partialBitPack,
                            unsigned char *bitPacked,
                            int size,
                            int atomSize,
                            int remainingSize,
                            int startingAtom) {
    if (remainingSize <= 0) {
        if (checkPresent(tree, partialBitPack, size) == NOTPRESENT) {
            /* Any section being not present is sufficient to determine non-presence. */
            return NOTPRESENT;
        }
        return PRESENT;
    }
    for (int i = startingAtom; i <= (tree->numPieces - remainingSize); i++) {
        writeNewBitsnCr(partialBitPack, atomSize * (size - remainingSize), bitPacked, atomSize * i, atomSize);
        int ppr = packPartialCheck(tree, partialBitPack, bitPacked, size, atomSize, remainingSize - 1, i + 1);
        if (ppr == NOTPRESENT) {
            return NOTPRESENT;
        }
    }
    /* Default behaviour - return present */
    return PRESENT;
}

/* Former inner: packPartial from insertRadixTreenCr */
static void packPartialInsert(struct radixTree *tree,
                              unsigned char *partialBitPack,
                              unsigned char *bitPacked,
                              int size,
                              int atomSize,
                              int remainingSize,
                              int startingAtom) {
    if (remainingSize <= 0) {
        if (checkPresent(tree, partialBitPack, size) == NOTPRESENT) {
            /* New value. */
            insertRadixTree(tree, partialBitPack, size);
        }
        return;
    }
    for (int i = startingAtom; i <= (tree->numPieces - remainingSize); i++) {
        writeNewBitsnCr(partialBitPack, atomSize * (size - remainingSize), bitPacked, atomSize * i, atomSize);
        packPartialInsert(tree, partialBitPack, bitPacked, size, atomSize, remainingSize - 1, i + 1);
    }
}

/* Checks if all state sections of length s are in the radix tree. */
int checkPresentnCr(struct radixTree *tree, unsigned char *bitPacked, int size) {
    int pBits = calcBits(tree->numPieces);
    int hBits = calcBits(tree->height);
    int wBits = calcBits(tree->width);
    int atomSize = pBits + hBits + wBits;

    /* Size * atoms, so bits contain location of size pieces. */
    int bitCount = atomSize * size;

    unsigned char *partialBitPack = (unsigned char *) calloc((bitCount + (BITS_PER_BYTE - 1) / BITS_PER_BYTE), sizeof(unsigned char));
    assert(partialBitPack);

    /* Start from each location */
    for(int i = 0; i <= (tree->numPieces - size); i++){
        writeNewBitsnCr(partialBitPack, 0, bitPacked, atomSize * i, atomSize);
        int ppr = packPartialCheck(tree, partialBitPack, bitPacked, size, atomSize, size - 1, i + 1);
        if(ppr == NOTPRESENT) {
            free(partialBitPack);
            return NOTPRESENT;
        }
    }
    assert((tree->numPieces - size) >= 0);
    free(partialBitPack);
    /* No missing atom combinations found. */
    return PRESENT;
}

/* Inserts sections of appropriate length */
void insertRadixTreenCr(struct radixTree *tree, unsigned char *bitPacked, int size) {
    int pBits = calcBits(tree->numPieces);
    int hBits = calcBits(tree->height);
    int wBits = calcBits(tree->width);
    int atomSize = pBits + hBits + wBits;

    /* Size * atoms, so bits contain location of size pieces. */
    int bitCount = atomSize * size;

    unsigned char *partialBitPack = (unsigned char *) calloc((bitCount + (BITS_PER_BYTE - 1) / BITS_PER_BYTE), sizeof(unsigned char));
    assert(partialBitPack);

    /* Start from each location */
    for(int i = 0; i <= (tree->numPieces - size); i++){
        writeNewBitsnCr(partialBitPack, 0, bitPacked, atomSize * i, atomSize);
        packPartialInsert(tree, partialBitPack, bitPacked, size, atomSize, size - 1, i + 1);
    }
    assert((tree->numPieces - size) >= 0);
    free(partialBitPack);
}

void freeRadixTree(struct radixTree *tree) {
    if(! tree) {
        return;
    }
    if(tree->prefixBytes) {
        free(tree->prefixBytes);
    }
    if(tree->prefixBitsBytes) {
        free(tree->prefixBitsBytes);
    }
    if(tree->prefixBitStartBytes) {
        free(tree->prefixBitStartBytes);
    }
    if(tree->branchABytes) {
        free(tree->branchABytes);
    }
    if(tree->branchBBytes) {
        free(tree->branchBBytes);
    }
    free(tree);
}
