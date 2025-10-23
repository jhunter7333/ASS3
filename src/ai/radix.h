/* 
 * Library for radix tree and related functions, bit-packed approach
 * which stores data in a space-compact method.
 * Written by Grady Fitzpatrick for COMP20003 Assignment 3 2025
*/
#ifndef __RADIX__
#define __RADIX__

#define NOTPRESENT (0)
#define PRESENT (1)

/* Index to use where no child is present - used to generate memory errors for access errors. */
#define NOCHILD (-1)

#include <stdint.h>
#include <unistd.h>

struct radixTree;

/* Return memory used in radix tree. */
int queryRadixMemoryUsage(struct radixTree *tree);

/* Helper utility to calculate the number of bits required to store a number */
int calcBits(int x);

/* Helper function. Gets the bit at bitIndex from the string s. */
int getBit(unsigned char *s, unsigned int bitIndex);

// set the i-th bit to true
void bitOn( unsigned char A[], unsigned int bitIndex );

// set the i-th bit to false
void bitOff( unsigned char A[], unsigned int bitIndex );

/*
	Creates new radix tree, numPieces is the number of pieces that are on the 
	board, height and width are of the board.
*/
struct radixTree *getNewRadixTree(int numPieces, int height, int width);

/* 
	Checks if the state is present in the radix tree. AtomCount is the number 
	of pieces in the bitPacked data.
*/
int checkPresent(struct radixTree *tree, unsigned char *bitPacked, int atomCount);

/* Inserts the state into the radix tree. */
void insertRadixTree(struct radixTree *tree, unsigned char *bitPacked, int atomCount);

/* Checks if all state sections of length s are in the radix tree. */
int checkPresentnCr(struct radixTree *tree, unsigned char *bitPacked, int size);

/* Inserts all sections of appropriate length s */
void insertRadixTreenCr(struct radixTree *tree, unsigned char *bitPacked, int size);

/* Free radix tree */
void freeRadixTree(struct radixTree *tree);

#endif
