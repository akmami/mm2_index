TARGET := segram
SRCS := $(wildcard *.c)
OBJS := $(SRCS:.c=.o)
CURRENT_DIR := $(shell pwd)

CC := cc
CFLAGS := -g -Wall -O2 -Wc++-compat -DHAVE_KALLOC
LDFLAGS := -lm -pthread -lz

# ========================================
#  Build Rules
# ========================================

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(CXXLIBS) $(LDFLAGS)
	rm -f $(OBJS)

%.o: %.c
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

# ========================================
#  Utility Targets
# ========================================

clean: 
	@echo "Cleaning"
	rm -f $(OBJS)
	rm -f $(TARGET)

install: clean $(TARGET)
