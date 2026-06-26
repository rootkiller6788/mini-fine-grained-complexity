CC = gcc
CFLAGS = -Wall -Wextra -O2 -I./include -std=c99
AR = ar

SRC = $(wildcard src/*.c)
OBJ = $(patsubst src/%.c,build/%.o,$(SRC))
LIB = build/libtsapsp.a

all: $(LIB) tests/run_tests

$(LIB): $(OBJ)
	@mkdir -p build
	$(AR) rcs $@ $^

build/%.o: src/%.c
	@mkdir -p build
	$(CC) $(CFLAGS) -c $< -o $@

tests/run_tests: tests/test_main.c $(LIB)
	$(CC) $(CFLAGS) -o $@ $< -Lbuild -ltsapsp -lm

test: tests/run_tests
	cd tests && ./run_tests

clean:
	rm -rf build/ tests/run_tests

.PHONY: all test clean
