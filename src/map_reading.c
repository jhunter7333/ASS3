/*
** EPITECH PROJECT, 2017
** PSU_my_sokoban_2017
** File description:
** Main function for the my_sokoban
** Edited by Thomas Minuzzo 2024 - adapted to Chessformer
** Edited by Grady Fitzpatrick 2025 - adapted to Impassable Gate
*/

#include <ncurses.h>
#include <unistd.h>
#include <stdlib.h>
#include <fcntl.h>
#include "../include/libmy.h"
#include "../include/gate.h"

#define READ_BUFFER_SIZE (10000)

char *open_map(char const *path) {
	int reading;
	char *buffer;

	reading = open(path, O_RDONLY);
	if (reading == -1) {
		write(2, "No such file or directory\n", 26);
		exit (84);
	}
	buffer = read_map(reading);
	close(reading);
	return (buffer);
}

char *read_map(int reading) {
	char *buffer = malloc(sizeof(char) * READ_BUFFER_SIZE);
	int size = 32;

	size = read(reading, buffer, READ_BUFFER_SIZE);
	if (size == -1) {
		exit(84);
	}
	buffer[size] = '\0';
	return (buffer);
}

int count_columns(gate_t gate, int position) {
	int columns = 0;

	for (; gate.buffer[position] != '\n'; position++) {
		columns++;
	}
	return (columns);
}

gate_t count_lines(gate_t gate) {
	gate.lines = 0;

	for (int i = 0; gate.buffer[i] != '\0'; i++) {
		if (gate.buffer[i] == '\n' || gate.buffer[i] == '\0')
			gate.lines++;
	}
	return (gate);
}

#define MAX_COLUMNS (26+2)
#define MAX_ROWS (9+2)
#include <assert.h>

gate_t make_map(char const *path, gate_t gate) {
	gate.buffer = open_map(path);
	gate.num_pieces = 0;
	gate = count_lines(gate);
	int k = 0;
	int columns = 0;
	gate.num_chars_map = 0;
	gate.map = malloc(sizeof(char *) * gate.lines);
	gate.map_save = malloc(sizeof (char *) * gate.lines);
	for (int j = 0; j < gate.lines; j++) {
		columns = count_columns(gate, k);
		gate.num_chars_map += columns;
		gate.map[j] = malloc(sizeof(char) * columns + 1);
		gate.map_save[j] = malloc(sizeof(char) * columns + 1);
		for (int i = 0; i < columns; i++) {
			gate.map[j][i] = gate.buffer[k];
			gate.map_save[j][i] = gate.buffer[k];
			gate.map[j][i+1] = '\0';
			gate.map_save[j][i+1] = '\0';
			k++;
		}
		k++;
	}
	assert(columns <= MAX_COLUMNS);
	assert(gate.lines <= MAX_ROWS);
	return (gate);
}
