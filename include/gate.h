/*
** EPITECH PROJECT, 2017
** bsq.h
** File description:
** Contain all the prototypes needed for BSQ
** Modified for COMP20003 Assignment 3 2025
*/

#ifndef BSQ_H
#define BSQ_H
	#define MAX_PIECES 10
	typedef struct gate {
		char *buffer; // Buffer for reading in the puzzle
		char **map; //A line by line map of chars representing the game state
		char **map_save; // A line-by-line map of chars, used for temporarily holding the state
		int lines; //The number of rows
		int player_x; //The player's current x location
		int player_y; //The player's current y location
		char const *base_path; // The file path of the puzzle
		char* soln; //To store the current moves made in this state
		int num_chars_map; //The number of chars total (divide by lines to get columns)
		int num_pieces; // The number of pieces on the board
		int piece_x[MAX_PIECES]; // x locations of part of each piece with 
								 // lowest y (tie-breaking with lowest x)
		int piece_y[MAX_PIECES]; // y locations of part of each piece with 
								 // lowest y (tie-breaking with lowest x)
	} gate_t;
	int helper(void);
	char *read_map(int reading);
	char *open_map(char const *path);
	gate_t make_map(char const *path, gate_t gate);
	int play(char const *path);
	gate_t count_lines(gate_t gate);
	int count_columns(gate_t gate, int position);
	gate_t check_if_player(gate_t gate, int y, int x);
	gate_t check_if_piece(gate_t gate, int y, int x, int piece);
	gate_t find_player(gate_t gate);
	gate_t find_pieces(gate_t gate);
	gate_t key_check(gate_t gate, char pieceNumber, char direction);
	gate_t attempt_move(gate_t gate, char pieceNumber, char direction);
	gate_t move_location(gate_t gate, char piece, char direction);
	int part_can_move(gate_t gate, int y, int x, char direction);
	void win_check(gate_t gate);
	void map_check(gate_t gate);
	int count_case_number(int y, int x, gate_t gate);
	int count_goal_square(int y, int x, gate_t gate);
	int count_player(int y, int x, gate_t gate);
	gate_t game_management(gate_t gate);
	int check_tile(int y, int x, gate_t gate);
#endif
