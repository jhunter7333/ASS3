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
#include <fcntl.h>
#include <stdlib.h>
#include "../include/libmy.h"
#include "../include/gate.h"
#include "ai/ai.h"

int main(int argc, char const **argv) {
	if (argc < 2 || argc > 3){
		helper();
		return (84);
	}
	if (argv[1][0] == '-' && argv[1][1] == 'h') {
		return(helper());
	} else if (argv[1][0] == '-' && argv[1][1] == 's') {
		solve(argv[2]);
		return 0;
	} else if (argv[1][0] != '-') {
		helper();
		return(play(argv[1]));
	}
	helper();
	return (84);
}
