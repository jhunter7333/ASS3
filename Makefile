##
## EPITECH PROJECT, 2017
## Makefile
## File description:
## Makefile
##

CC	=	gcc -Wall -Wextra -O3 -g


RM	=	rm -f

NAME	=	gate

SRC	=	src/main.c	\
		src/helper.c	\
		src/key_check.c	\
		src/find_player.c	\
		src/map_check.c	\
		src/map_reading.c	\
		src/movement.c	\
		src/play.c	\
		src/win_check.c	\
		lib/my_putchar.c	\
		lib/my_putstr.c	\
		src/ai/radix.o \
		src/ai/ai.o \
		src/ai/utils.o

CFLAGS	+=	-I./include/

OBJ	=	$(SRC:.c=.o)

all:	$(NAME)

$(NAME):	$(OBJ)
	$(CC) -o $(NAME) $(OBJ) -lncurses

clean:
	$(RM) $(OBJ)

fclean: clean
	$(RM) $(NAME)

re:	fclean all

runmanual:
	make
	./gate test_puzzles/capability1
	./gate test_puzzles/capability2
	./gate test_puzzles/capability3
	./gate test_puzzles/capability4
	./gate test_puzzles/capability5
	./gate test_puzzles/capability6
	./gate test_puzzles/capability7
	./gate test_puzzles/capability8
	./gate test_puzzles/capability9
	./gate test_puzzles/capability10
	./gate test_puzzles/capability11
	./gate test_puzzles/capability12
	./gate test_puzzles/capability13
	./gate test_puzzles/impassable1
	./gate test_puzzles/impassable2
	./gate test_puzzles/impassable3

runtests:
	make clean
	make
	./gate -s test_puzzles/capability1
	./gate -s test_puzzles/capability2
	./gate -s test_puzzles/capability3
	./gate -s test_puzzles/capability4
	./gate -s test_puzzles/capability5
	./gate -s test_puzzles/capability6
	./gate -s test_puzzles/capability7
	./gate -s test_puzzles/capability8
	./gate -s test_puzzles/capability9
	./gate -s test_puzzles/capability10
	./gate -s test_puzzles/capability11
	./gate -s test_puzzles/capability12
	./gate -s test_puzzles/capability13
	./gate -s test_puzzles/impassable1
	./gate -s test_puzzles/impassable2
	./gate -s test_puzzles/impassable3

.PHONY: all clean fclean re
