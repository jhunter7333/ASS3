/*
** EPITECH PROJECT, 2017
** PSU_my_sokoban_2017
** File description:
** Function that manage key press for sokoban
** Edited by Thomas Minuzzo 2024 - adapted to Chessformer
** Edited by Grady Fitzpatrick 2025 - adapted to Impassable Gate
*/

#include <ncurses.h>
#include <unistd.h>
#include <stdlib.h>
#include <fcntl.h>
#include "../include/libmy.h"
#include "../include/gate.h"

gate_t key_check(gate_t gate, char pieceNumber, char direction) {
	gate = move_location(gate, pieceNumber, direction);
	return (gate);
}

gate_t attempt_move(gate_t gate, char pieceNumber, char direction) {
	gate = move_location(gate, pieceNumber, direction);
	return (gate);
}
