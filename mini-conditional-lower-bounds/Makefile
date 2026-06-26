CC       = gcc
CFLAGS   = -std=c99 -Wall -Wextra -O2 -g
LDFLAGS  = -lm
INC_DIR  = include
SRC_DIR  = src
TEST_DIR = tests
BUILD_DIR = build

SRCS = $(SRC_DIR)/condlb.c \
       $(SRC_DIR)/seth_bound.c \
       $(SRC_DIR)/ov_reduction.c \
       $(SRC_DIR)/string_hardness.c \
       $(SRC_DIR)/threesum_apsp.c \
       $(SRC_DIR)/reduction_web.c

OBJS = $(patsubst $(SRC_DIR)/%.c,$(BUILD_DIR)/%.o,$(SRCS))
TEST_SRC = $(TEST_DIR)/test_main.c
TEST_EXE = $(BUILD_DIR)/test_condlb

all: $(TEST_EXE)

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.c | $(BUILD_DIR)
	$(CC) $(CFLAGS) -I$(INC_DIR) -c $< -o $@

$(TEST_EXE): $(OBJS) $(TEST_SRC)
	$(CC) $(CFLAGS) -I$(INC_DIR) $(TEST_SRC) $(OBJS) $(LDFLAGS) -o $@

test: $(TEST_EXE)
	./$(TEST_EXE)

clean:
	rm -rf $(BUILD_DIR)

.PHONY: all test clean
