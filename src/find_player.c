/*
** EPITECH PROJECT, 2017
** PSU_my_sokoban_2017
** File description:
** Function for finding the player on the map
** Edited by Thomas Minuzzo 2024 - adapted to Chessformer
** Edited by Grady Fitzpatrick 2025 - adapted to Impassable Gate
*/

#include <ncurses.h>
#include <unistd.h>
#include <stdlib.h>
#include <fcntl.h>
#include "../include/libmy.h"
#include "../include/gate.h"

gate_t find_player(gate_t gate) {
	gate.player_x = 0;
	gate.player_y = 0;
	for (int i = 0; i < gate.lines; i++) {
		for (int j = 0; gate.map[i][j] != '\0'; j++) {
			gate = check_if_player(gate, i, j);
		}
	}
	return (gate);
}

gate_t find_pieces(gate_t gate) {
	for(int i = 0; i < MAX_PIECES; i++) {
		gate.piece_x[i] = -1;
		gate.piece_y[i] = -1;
	}
	for (int i = 0; i < gate.lines; i++) {
		for (int j = 0; gate.map[i][j] != '\0'; j++) {
			if(! ((gate.map[i][j] >= '0' && gate.map[i][j] <= '9') || (gate.map[i][j] >= 'H' && gate.map[i][j] <= 'Q'))){
				// Not a piece.
				continue;
			}
			int piece = gate.map[i][j];
			int letterPiece = gate.map[i][j];
			if(piece >= '0' && piece <= '9') {
				letterPiece = piece - '0' + 'H';
			} else {
				piece = letterPiece - 'H' + '0';
			}
			int pieceIndex = piece - '0';
			if((pieceIndex + 1) > gate.num_pieces) {
				gate.num_pieces = pieceIndex + 1;
			}
			// Only update if we find a piece for the first time.
			if(gate.piece_x[pieceIndex] == -1) {
				gate = check_if_piece(gate, i, j, piece);
			}
		}
	}
	return (gate);
}

gate_t check_if_player(gate_t gate, int y, int x) {
	if (gate.map[y][x] == '0' || gate.map[y][x] == 'H') {
		gate.player_x = x;
		gate.player_y = y;
	}
	return (gate);
}

gate_t check_if_piece(gate_t gate, int y, int x, int piece) {
	if (gate.map[y][x] == (piece) || gate.map[y][x] == ('H' + piece - '0')) {
		gate.piece_x[piece - '0'] = x;
		gate.piece_y[piece - '0'] = y;
	}
	return (gate);
}

