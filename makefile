CC=gcc
CFLAGS=-c
LFLAGS=
TARGET=comp_demo

SRCS=$(wildcard *.c)
OBJS=$(patsubst %.c,%.o,$(SRCS))

all: $(TARGET)

debug: CFLAGS += -g -Wall
debug: $(TARGET)

%.o: %.c
	$(CC) $(CFLAGS) $< -o $@

$(TARGET): $(OBJS)
	$(CC) $(LFLAGS) $^ -o $@

clean:
	rm *.o
	rm $(TARGET)
