#include <time.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdbool.h>
#include <assert.h>

#include "ai.h"
#include "gate.h"
#include "radix.h"
#include "utils.h"

#define DEBUG 0

#define UP 'u'
#define DOWN 'd'
#define LEFT 'l'
#define RIGHT 'r'

char directions[] = {UP, DOWN, LEFT, RIGHT};
char invertedDirections[] = {DOWN, UP, RIGHT, LEFT};
char pieceNames[] = {'0', '1', '2', '3', '4', '5', '6', '7', '8', '9'};

// Queue Structure
typedef struct QueueNode {
    gate_t *state;
    char   *path;
    int     depth;
    struct QueueNode *next;
} QueueNode;

static QueueNode *q_head = NULL;
static QueueNode *q_tail = NULL;

// Append Queue 
static void queue_push(QueueNode *n) { 
	n->next = NULL; 
	if(!q_tail){ 
		q_head = q_tail = n; 
	}else { 
		q_tail->next = n; 
		q_tail = n; 
	} 
}

// Pop out top of Queue 
static QueueNode* q_pop(void) { 
	if(!q_head){
		return NULL; 
	}

	QueueNode *n=q_head; 
	q_head=n->next; 

	if(!q_head) {
		q_tail=NULL; return n; 
	}

	return n;
}



/* Compare two maps to detect if move changed positions/map */
int maps_equal(gate_t *a, gate_t *b) {
	if (a->lines != b->lines){
		return 0;
	}

	for (int i = 0; i < a->lines; i++) {
		if (strcmp(a->map[i], b->map[i]) != 0){
			return 0;
		}
	}

	return 1;
}

/* Append one step to existing path */
static char* append_path(const char *parent, char piece, char dir) {
	size_t parent_length = 0;

	if (parent){
		parent_length = strlen(parent);
	}

	char *result = (char*)malloc(parent_length + 2 + 1);
	assert(result);

	if (parent){
		memcpy(result, parent, parent_length);
	}

	result[parent_length] = piece;
	result[parent_length + 1] = dir;
	result[parent_length + 2] = '\0';
	return result;
	}

/**
 * Given a game state, work out the number of bytes required to store the state.
*/
int getPackedSize(gate_t *gate);

/**
 * Store state of puzzle in map.
*/
void packMap(gate_t *gate, unsigned char *packedMap);

/**
 * Check if the given state is in a won state.
 */
bool winning_state(gate_t gate);

/* duplicate_state
       - Allocate a new gate_t and copy the state.
       - Copies all scalar/array fields through the struct.
       - Copies the 2D maps row-by-row.
       - Returns a pointer to the duplicated state or NULL if fails.
*/
gate_t* duplicate_state(gate_t* gate) {
	gate_t* duplicate = (gate_t*)malloc(sizeof(gate_t));
	assert(duplicate);

	// Copies scalars and pointers
	*duplicate = *gate;
	duplicate->buffer = NULL;

	if (gate->lines > 0) {
        duplicate->map = (char**)malloc(sizeof(char*) * gate->lines);

        duplicate->map_save = (char**)malloc(sizeof(char*) * gate->lines);
		assert(duplicate->map && duplicate->map_save);

        /* Copy each row of map and map_save */
        for (int i = 0; i < gate->lines; i++) {
            /* map rows */
            if (gate->map && gate->map[i]) {
                size_t n = strlen(gate->map[i]) + 1;
                duplicate->map[i] = (char*)malloc(n);
                memcpy(duplicate->map[i], gate->map[i], n);
            } else {
                duplicate->map[i] = NULL;
            }

            /* map_save row */
            if (gate->map_save && gate->map_save[i]) {
                size_t m = strlen(gate->map_save[i]) + 1;
                duplicate->map_save[i] = (char*)malloc(m);
                memcpy(duplicate->map_save[i], gate->map_save[i], m);
            } else {
                duplicate->map_save[i] = NULL;
            }
        }
    } else {
        duplicate->map = NULL;
        duplicate->map_save = NULL;
    }

	return duplicate;

}

/* free_state
    - Release all dynamic memory in search state.
    - Frees each row of map and map_save,  pointer arrays, and the state structure.
*/
void free_state(gate_t* stateToFree, gate_t *init_data) {
	
	if (!stateToFree) return;

    /* free all map rows */
    if (stateToFree->map) {
        for (int i = 0; i < stateToFree->lines; i++) {
            free(stateToFree->map[i]);
            stateToFree->map[i] = NULL;
        }
        free(stateToFree->map);
        stateToFree->map = NULL;
    }

	if (stateToFree->map_save) {
        for (int i = 0; i < stateToFree->lines; i++) {
			free(stateToFree->map_save[i]);
            stateToFree->map_save[i] = NULL;	
        }
		free(stateToFree->map_save);
        stateToFree->map_save = NULL;
    }

    free(stateToFree);

}


/* free_initial_state
    - Releases the dynamic memory owned by the initial state.
    - Includes: buffer, map rows, map_save rows and pointer arrays.
*/
void free_initial_state(gate_t *init_data) {
	
	/* Free map and map_save rows and arrays */
    for (int i = 0; i < init_data->lines; i++) {
        free(init_data->map[i]);
		free(init_data->map_save[i]);
        init_data->map[i] = NULL;
		init_data->map_save[i] = NULL;
    }
	
    free(init_data->map);
    free(init_data->map_save);
	init_data->map = NULL;
    init_data->map_save = NULL;

	// Free buffer
	free(init_data->buffer);
    init_data->buffer = NULL;

}


static QueueNode* apply_action(const QueueNode *parent, char piece, char dir, gate_t *init_data) {
    gate_t *child_state = duplicate_state(parent->state);
    *child_state = attempt_move(*child_state, piece, dir);

    if (maps_equal(parent->state, child_state)) {
        free_state(child_state, init_data);
        return NULL;
    }

    QueueNode *child = (QueueNode*)malloc(sizeof *child);
    assert(child);
    child->state = child_state;
    child->path  = append_path(parent->path, piece, dir);
    child->depth = parent->depth + 1;
    child->next  = NULL;
    return child;
}


/**
 * Find a solution by exploring all possible paths
 */
void find_solution(gate_t* init_data) {
	/* Location for packedMap. */
	int packedBits = getPackedSize(init_data);
	int packedBytes = (packedBits+ 7) / 8;
	unsigned char *packedMap = (unsigned char *) calloc(packedBytes, sizeof(unsigned char));
	assert(packedMap);

	bool has_won = false;
	int dequeued = 0;
	int enqueued = 0;
	int duplicatedNodes = 0;
	char *soln = "";
	double start = now();
	double elapsed;
	
	// Algorithm 1 is a width n + 1 search
	int w = init_data->num_pieces + 1;

	/* ----------- ALGORITHMS ----------- */

	int mapHeight = init_data->lines;
	int mapWidth = (init_data->num_chars_map / init_data->lines);
	int mapPieces = init_data->num_pieces;

	const int pBits = calcBits(init_data->num_pieces);
    const int hBits = calcBits(init_data->lines);
    const int wBits = calcBits(mapWidth);

	/* ALG 2
	struct radixTree *radixTree = getNewRadixTree(mapPieces, mapHeight, mapWidth);
	assert(radixTree);
	*/ /* Comment out for ALG 1 and 3 */

	/* ALG3 */
	struct radixTree **rts = (struct radixTree**)calloc(w + 1, sizeof *rts);
	assert(rts);
	/* Comment out for ALG 1 and 2 */

	/* ALG 3 FOR LOOP */
	for (int k = 1; k <= w && !has_won; k++){
	/* ALG 3 FOR LOOP */

		for (int s = 1; s <= k && s <= mapPieces; s++){
			if (rts[s]){
				freeRadixTree(rts[s]);
				rts[s] = NULL;
			} 
			rts[s] = getNewRadixTree(mapPieces, mapHeight, mapWidth);
			assert(rts[s]);
		}

		q_head = NULL; 
		q_tail = NULL;

		/* Initialise Root */
		QueueNode *root = (QueueNode*)malloc(sizeof *root);
		root->state = duplicate_state(init_data);
		root->path  = strdup("");
		root->depth = 0;
		root->next  = NULL;

		/* Enqueue root */
		queue_push(root);
		enqueued++;

		packMap(root->state, packedMap);
		for (int s = 1; s <= k && s <= mapPieces; s++) {
			insertRadixTreenCr(rts[s], packedMap, s);
		}

		/* ALG 2
		packMap(root->state, packedMap);
		insertRadixTree(radixTree, packedMap, init_data->num_pieces);
		*/ /* ALG 2 Comment out for ALG 1 and 3 */


		/* Until queue is empty */
		while (q_head) {

			QueueNode *n = q_pop();
			dequeued++;

			if (winning_state(*n->state)) {

				soln = strdup(n->path);  /* Save a copy to print after loop */
				free_state(n->state, init_data);
				free(n->path);
				free(n);

				/* Drain remaining queue nodes */
				while (q_head) {
					QueueNode *t = q_pop();
					free_state(t->state, init_data);
					free(t->path);
					free(t);
				}

				has_won = true;
				w = k; // ALG 3 ONLY
				break;
			}

			/* for each move action */
			for (int p = 0; p < init_data->num_pieces; p++) {
				char piece = pieceNames[p];

				for (int i = 0; i < 4; i++) {
					char dir = directions[i];

					/* Apply action and create new child */
					QueueNode *child = apply_action(n, piece, dir, init_data);
					if (!child) continue;

					/* ALG 2 
					packMap(child->state, packedMap);

					if (checkPresent(radixTree, packedMap, mapPieces) == PRESENT) {
						duplicatedNodes++;	
						free_state(child->state, init_data);
						free(child->path);
						free(child);
						continue;
					}

					insertRadixTree(radixTree, packedMap, mapPieces);
					*/ 
					/* ALG 2 Comment out for ALG 1 and 3 */


					/* ALG 3 ONLY */
					int novel = 0;
					packMap(child->state, packedMap);

					 /* for each size, if section not in radixTrees[s], novel is true */
					for (int s = 1; s <= k && s <= mapPieces; s++) {
						if (checkPresentnCr(rts[s], packedMap, s) == NOTPRESENT) {
							novel = 1;
							break;
                    	}
                	}

					/* if not novel then duplicatedNode */
					if (!novel) {
						duplicatedNodes++;
						free_state(child->state, init_data);
						free(child->path);
    					free(child);
						continue;
					}

					for (int s = 1; s <= k && s <= mapPieces; s++) {
						insertRadixTreenCr(rts[s], packedMap, s);
					}
					/* ALG 3 Comment out for ALG 1 and 2 */

					enqueued++;    
					queue_push(child);
				}
			}

			/* finished with current node */
			free_state(n->state, init_data);
			free(n->path);
			free(n);
		}

	// ALG 3 ONLY 
	}
	// ALG 3 ONLY

	 /* Output statistics */
	elapsed = now() - start;
	printf("Solution path: ");
	printf("%s\n", soln);
	printf("Execution time: %lf\n", elapsed);
	printf("Expanded nodes: %d\n", dequeued);
	printf("Generated nodes: %d\n", enqueued);
	printf("Duplicated nodes: %d\n", duplicatedNodes);
	int memoryUsage = 0;

	/* ALG 2 ONLY
	memoryUsage += queryRadixMemoryUsage(radixTree);
	freeRadixTree(radixTree);
	*/ /* ALG 2 Comment out for ALG 1 and 3 */

	// Algorithm 3: Memory usage, uncomment to add.
	for(int i = 0; i <= w; i++) {
		if (rts[i]){
			memoryUsage += queryRadixMemoryUsage(rts[i]);
		} 
	}

	/* Free novelty structures */
	for (int s = 1; s <= mapPieces; s++) {
		if (rts[s]) {
			freeRadixTree(rts[s]);
		}
	}
	free(rts);

	/* ALG 3 Comment out for ALG 1 and 2 */


	printf("Auxiliary memory usage (bytes): %d\n", memoryUsage);
	printf("Number of pieces in the puzzle: %d\n", init_data->num_pieces);
	printf("Number of steps in solution: %ld\n", strlen(soln)/2);
	int emptySpaces = 0;

	for (int row = 0; row < init_data->lines; row++) {
		for (int column = 0; init_data->map[row][column] != '\0'; column++) {
			if (init_data->map[row][column] == ' ') {
				emptySpaces++;
			}
		}
	}
	
	printf("Number of empty spaces: %d\n", emptySpaces);
	printf("Solved by IW(%d)\n", w);
	printf("Number of nodes expanded per second: %lf\n", (dequeued + 1) / elapsed);

	/* Free associated memory. */
	if(packedMap) {
		free(packedMap);
	}
	/* Free initial map. */
	free_initial_state(init_data);
	free(soln);
}

/**
 * Given a game state, work out the number of bytes required to store the state.
*/
int getPackedSize(gate_t *gate) {
	int pBits = calcBits(gate->num_pieces);
    int hBits = calcBits(gate->lines);
    int wBits = calcBits(gate->num_chars_map / gate->lines);
    int atomSize = pBits + hBits + wBits;
	int bitCount = atomSize * gate->num_pieces;
	return bitCount;
}

/**
 * Store state of puzzle in map.
*/
void packMap(gate_t *gate, unsigned char *packedMap) {
	int pBits = calcBits(gate->num_pieces);
    int hBits = calcBits(gate->lines);
    int wBits = calcBits(gate->num_chars_map / gate->lines);
	int bitIdx = 0;
	for(int i = 0; i < gate->num_pieces; i++) {
		for(int j = 0; j < pBits; j++) {
			if(((i >> j) & 1) == 1) {
				bitOn( packedMap, bitIdx );
			} else {
				bitOff( packedMap, bitIdx );
			}
			bitIdx++;
		}
		for(int j = 0; j < hBits; j++) {
			if(((gate->piece_y[i] >> j) & 1) == 1) {
				bitOn( packedMap, bitIdx );
			} else {
				bitOff( packedMap, bitIdx );
			}
			bitIdx++;
		}
		for(int j = 0; j < wBits; j++) {
			if(((gate->piece_x[i] >> j) & 1) == 1) {
				bitOn( packedMap, bitIdx );
			} else {
				bitOff( packedMap, bitIdx );
			}
			bitIdx++;
		}
	}
}

/**
 * Check if the given state is in a won state.
 */
bool winning_state(gate_t gate) {
	for (int i = 0; i < gate.lines; i++) {
		for (int j = 0; gate.map_save[i][j] != '\0'; j++) {
			if (gate.map[i][j] == 'G' || (gate.map[i][j] >= 'I' && gate.map[i][j] <= 'Q')) {
				return false;
			}
		}
	}
	return true;
}

void solve(char const *path)
{
	/**
	 * Load Map
	*/
	gate_t gate = make_map(path, gate);
	
	/**
	 * Verify map is valid
	*/
	map_check(gate);

	/**
	 * Locate player x, y position
	*/
	gate = find_player(gate);

	/**
	 * Locate each piece.
	*/
	gate = find_pieces(gate);
	
	gate.base_path = path;

	find_solution(&gate);

}
