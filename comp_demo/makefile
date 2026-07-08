CC=gcc
CFLAGS=-c
LFLAGS=
TARGET=comp_demo

OBJDIR=obj

SRCS=$(wildcard *.c)
OBJS=$(patsubst %.c,$(OBJDIR)/%.o,$(SRCS))

all: $(TARGET)

debug: CFLAGS += -g -Wall
debug: $(TARGET)

$(OBJDIR)/%.o: %.c | obj
	$(CC) $(CFLAGS) $< -o $@

$(TARGET): $(OBJS)
	$(CC) $^ $(LFLAGS) -o $@

obj:
	mkdir -p $@
	mkdir -p outputs

clean:
	rm -r $(OBJDIR)
	rm $(TARGET)
