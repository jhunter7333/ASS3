/*
** EPITECH PROJECT, 2017
** PSU_my_sokoban_2017
** File description:
** Function used to check if a map is valid or not
** Edited by Thomas Minuzzo 2024 - adapted to Chessformer
** Edited by Grady Fitzpatrick 2025 - adapted to Impassable Gate
*/

#include <unistd.h>
#include <stdlib.h>
#include <fcntl.h>
#include <assert.h>
#include "../include/libmy.h"
#include "../include/gate.h"

void map_check(gate_t gate) {
	int player = 0;
	int goal_squares = 0;
	for (int i = 0; i < gate.lines; i++) {
		for (int j = 0; gate.map_save[i][j] != '\0'; j++) {
			check_tile(i, j, gate);
			player += count_player(i, j, gate);
			goal_squares += count_goal_square(i, j, gate);
		}
	}
	assert(player > 0);
	assert(goal_squares > 0);
	assert(player == goal_squares);
}

int check_tile(int y, int x, gate_t gate) {
	// Avaliable characters are:
	// space - empty space.
	// 0 - 9 - part of a block.
	// G - Q - part of a tile or goal location.
	// # - wall.
	if ((! (gate.map_save[y][x] >= '0' && gate.map_save[y][x] <= '9'))
		&& (! (gate.map_save[y][x] >= 'G' && gate.map_save[y][x] <= 'Q'))
		&& gate.map_save[y][x] != '#' && gate.map_save[y][x] != ' '
		&& gate.map_save[y][x] != '\n') {
		write(2, "Unknown read character in map\n", 26);
		exit(84);
	}
	return (0);
}

int count_player(int y, int x, gate_t gate) {
	int i = 0;

	if (gate.map_save[y][x] == 'H' || gate.map_save[y][x] == '0') {
		i++;
	}
	return (i);
}

int count_goal_square(int y, int x, gate_t gate) {
	int i = 0;

	if (gate.map_save[y][x] >= 'G' && gate.map_save[y][x] <= 'Q') {
		i++;
	}
	return (i);
}

