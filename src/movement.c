/*
** EPITECH PROJECT, 2017
** PSU_my_sokoban_2017
** File description:
** Function that make the player move
** Edited by Thomas Minuzzo 2024 - adapted to Chessformer
** Edited by Grady Fitzpatrick 2025 - adapted to Impassable Gate
*/


#include <ncurses.h>
#include <unistd.h>
#include <stdlib.h>
#include <fcntl.h>
#include "../include/libmy.h"
#include "../include/gate.h"

// Tries to apply any legal action
gate_t move_location(gate_t gate, char piece, char direction){
	//Overwrite print above...
	mvprintw(gate.lines + 4, 0,
		"                                                                          ");

	// Valid move:
	if(!(direction == 'u' || 
		 direction == 'd' || 
		 direction == 'l' || 
		 direction == 'r')
	  ) {
		mvprintw(gate.lines + 4, 0, 
			"Illegal movement direction: %c not one of {u, d, l, r}", direction);
		return gate;
	}

	// Valid piece:
	if(!(piece >= '0' || piece <= '9')) {
		mvprintw(gate.lines + 4, 0, 
			"Illegal piece selection: %c not one of \[0, ..., 9]", piece);
		return gate;
	}

	// Find piece:
	// Representation of piece on goal location.
	char letterPiece = piece - '0' + 'H';
	// Reset location so that we get the right corner piece.
	gate.piece_x[piece - '0'] = -1;
	gate.piece_y[piece - '0'] = -1;
	for (int i = 0; i < gate.lines; i++) {
		for (int j = 0; gate.map[i][j] != '\0'; j++) {
			gate = check_if_piece(gate, i, j, piece);
			if (gate.piece_x[piece - '0'] != -1) {
				break;
			}
		}
		if (gate.piece_x[piece - '0'] != -1) {
			break;
		}
	}
	if (gate.piece_x[piece - '0'] == -1) {
		// Didn't find piece.
		mvprintw(gate.lines + 4, 0, "Unable to find piece on the board: %c", piece);
		return gate;
	}

	// Dry run:
	// Initially assume can move and look for a violation of that assumption.
	int can_move = 1;
	for (int i = gate.piece_y[piece - '0']; i < gate.lines; i++) {
		for (int j = 0; gate.map[i][j] != '\0'; j++) {
			if ((gate.map[i][j] == piece || gate.map[i][j] == letterPiece) 
				&& (! part_can_move(gate, i, j, direction))){
				can_move = 0;
				break;
			}
		}
	}
	if (! can_move) {
		// Didn't find piece.
		mvprintw(gate.lines + 4, 0, "Unable to move piece %c in direction %c", piece, direction);
		return gate;
	}

	// Make move:
	mvprintw(gate.lines + 4, 0, "Making move: %c%c", piece, direction);
	// Make backup of map.
	for (int i = 0; i < gate.lines; i++) {
		for (int j = 0; gate.map[i][j] != '\0'; j++) {
			gate.map_save[i][j] = gate.map[i][j];
		}
	}
	for (int i = 0; i < gate.lines; i++) {
		for (int j = 0; gate.map[i][j] != '\0'; j++) {
			if(gate.map_save[i][j] == piece || gate.map_save[i][j] == letterPiece) {
				if (direction == 'u'){
					if (gate.map_save[i - 1][j] == 'G' || gate.map_save[i - 1][j] == letterPiece){
						gate.map[i - 1][j] = letterPiece;
					} else {
						gate.map[i - 1][j] = piece;
					}
				} else if (direction == 'd') {
					if (gate.map_save[i + 1][j] == 'G' || gate.map_save[i + 1][j] == letterPiece){
						gate.map[i + 1][j] = letterPiece;
					} else {
						gate.map[i + 1][j] = piece;
					}
				} else if (direction == 'l') {
					if (gate.map_save[i][j - 1] == 'G' || gate.map_save[i][j - 1] == letterPiece){
						gate.map[i][j - 1] = letterPiece;
					} else {
						gate.map[i][j - 1] = piece;
					}
				} else if (direction == 'r') {
					if (gate.map_save[i][j + 1] == 'G' || gate.map_save[i][j + 1] == letterPiece){
						gate.map[i][j + 1] = letterPiece;
					} else {
						gate.map[i][j + 1] = piece;
					}
				}
				// Resolve piece on goal:
				// 1. Other piece moves onto location, no change.
				// 2. Otherwise becomes empty goal location.
				if (gate.map_save[i][j] == letterPiece && 
					(! (direction == 'u' && (gate.map_save[i + 1][j] == piece || gate.map_save[i + 1][j] == letterPiece))) &&
					(! (direction == 'd' && (gate.map_save[i - 1][j] == piece || gate.map_save[i - 1][j] == letterPiece))) &&
					(! (direction == 'l' && (gate.map_save[i][j + 1] == piece || gate.map_save[i][j + 1] == letterPiece))) &&
					(! (direction == 'r' && (gate.map_save[i][j - 1] == piece || gate.map_save[i][j - 1] == letterPiece)))
				   ) {
					gate.map[i][j] = 'G';
				}
				// Resolve piece on blank space:
				// Same logic.
				if (gate.map_save[i][j] == piece && 
					(! (direction == 'u' && (gate.map_save[i + 1][j] == piece || gate.map_save[i + 1][j] == letterPiece))) &&
					(! (direction == 'd' && (gate.map_save[i - 1][j] == piece || gate.map_save[i - 1][j] == letterPiece))) &&
					(! (direction == 'l' && (gate.map_save[i][j + 1] == piece || gate.map_save[i][j + 1] == letterPiece))) &&
					(! (direction == 'r' && (gate.map_save[i][j - 1] == piece || gate.map_save[i][j - 1] == letterPiece)))
				   ) {
					gate.map[i][j] = ' ';
				}
			}
		}
	}

	// Update piece.
	gate.piece_x[piece - '0'] = -1;
	gate.piece_y[piece - '0'] = -1;
	for (int i = 0; i < gate.lines; i++) {
		for (int j = 0; gate.map[i][j] != '\0'; j++) {
			gate = check_if_piece(gate, i, j, piece);
			if (gate.piece_x[piece - '0'] != -1) {
				break;
			}
		}
		if (gate.piece_x[piece - '0'] != -1) {
			break;
		}
	}

	return gate;
}

int part_can_move(gate_t gate, int y, int x, char direction) {
	int piece = gate.map[y][x];
	int letterPiece = gate.map[y][x];
	if(piece >= '0' && piece <= '9') {
		letterPiece = piece - '0' + 'H';
	} else {
		piece = letterPiece - 'H' + '0';
	}
	int can_move = 0;
	if (direction == 'u') {
		if(gate.map[y - 1][x] == ' ' 
			|| gate.map[y - 1][x] == 'G' || gate.map[y - 1][x] == letterPiece 
			|| gate.map[y - 1][x] == piece) {
			can_move = 1;
		}
	} else if (direction == 'd') {
		if(gate.map[y + 1][x] == ' ' 
			|| gate.map[y + 1][x] == 'G' || gate.map[y + 1][x] == letterPiece 
			|| gate.map[y + 1][x] == piece) {
			can_move = 1;
		}
	} else if (direction == 'l') {
		if(gate.map[y][x - 1] == ' ' 
			|| gate.map[y][x - 1] == 'G' || gate.map[y][x - 1] == letterPiece 
			|| gate.map[y][x - 1] == piece) {
			can_move = 1;
		}
	} else if (direction == 'r') {
		if(gate.map[y][x + 1] == ' ' 
			|| gate.map[y][x + 1] == 'G' || gate.map[y][x + 1] == letterPiece 
			|| gate.map[y][x + 1] == piece) {
			can_move = 1;
		}
	}
	return can_move;
}

