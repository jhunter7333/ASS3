/*
** EPITECH PROJECT, 2017
** PSU_my_sokoban_2017
** File description:
** Main function for the my_sokoban
** Edited by Thomas Minuzzo 2024 - adapted to Chessformer
** Edited by Grady Fitzpatrick 2025 - adapted to Impassable Gate
*/

#include <unistd.h>
#include "../include/libmy.h"

int helper(void) {
	my_putstr("USAGE\n");
	my_putstr("	./gate <-s> puzzle\n\n");
	my_putstr("DESCRIPTION\n");
	my_putstr(" Arguments within <> are optional\n");
	my_putstr("    -s                 calls the AI solver\n");
	return (0);
}
