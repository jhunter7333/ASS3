/*
** EPITECH PROJECT, 2017
** PSU_my_sokoban_2017
** File description:
** Function that manage the game
** ----------
** Adapted by Nir Lipo, 2021
** Edited by Thomas Minuzzo 2024 - adapted to Chessformer
** Edited by Grady Fitzpatrick 2025 - adapted to Impassable Gate
*/


#include <ncurses.h>
#include <unistd.h>
#include <stdlib.h>
#include <fcntl.h>
#include "../include/libmy.h"
#include "../include/gate.h"

/*********
* MACROS *
*********/
#include <string.h>
#define TERMINAL_TYPE (strcmp(getenv("TERM"), "xterm") == 0 ? "rxvt" : \
  getenv("TERM"))
SCREEN *mainScreen = NULL;
WINDOW *mainWindow = NULL;

int play(char const *path) {
	/**
	 * Load Map
	*/
	gate_t gate = make_map(path, gate);
	
	/**
	 * Count number of pieces and piece locations
	 * to verify the map is valid.
	*/
	map_check(gate);

	/**
	 * Locate player x, y position
	*/
	gate = find_player(gate);
	
	gate.base_path = path;

	mainScreen = newterm(TERMINAL_TYPE, stdout, stdin);
	set_term(mainScreen);
	int cols = 1;
	for(int i = 0; i < gate.lines; i++){
		if(strlen(gate.map[i]) > (size_t) cols){
			cols = strlen(gate.map[i]);
		}
	}
	mainWindow = newwin(gate.lines, cols, 0, 0);

	cbreak();
	noecho();
	keypad(stdscr, TRUE);
	clear();
	
	while (1) {
			refresh();
			// Print board
			int columns = gate.num_chars_map / gate.lines;
			char letter[2] = "0";
			for(int i = 1; i < columns - 1; i++){
				mvprintw(gate.lines, i, letter);
				letter[0] ++; //Move letter along, by editing string
			}
			char number[2] = "0";
			for (int i = 0; i < gate.lines; i++){
				if(i != 0 && i != gate.lines - 1){
					mvprintw(i, columns, number); //-3 for new line and null byte and array zero offset
					//mvprintw(i, gate.lines-3, number); //-3 for new line and null byte and array zero offset
					number[0]++; //Move number along, by editing string
				}
				mvprintw(i, 0, gate.map[i]);
			}
			move(gate.player_y, gate.player_x);

			usleep(500000);
			win_check(gate);// extra win check for gravity
			//Check for key presses
			gate = game_management(gate);	
	}
}

gate_t game_management(gate_t gate) {
	int piece = 0 , direction = 0;

	mvprintw(gate.lines + 1, 0, "Enter a piece number and a letter direction to move that piece in that direction");
	mvprintw(gate.lines + 2, 0, "Piece:     ");
	mvprintw(gate.lines + 3, 0, "Direction: ");

	piece = getch();
	mvprintw(gate.lines + 2, 0, "Piece: %c", piece);
	direction = getch();
	mvprintw(gate.lines + 3, 0, "Direction: %c", direction);

	gate = key_check(gate, piece, direction);
	win_check(gate);
	return (gate);
}
