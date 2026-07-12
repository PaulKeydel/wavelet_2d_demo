CC=gcc
CFLAGS=-c
LFLAGS=
TARGET1=comp_demo
TARGET2=reco_demo

OBJDIR=obj
OUTDIR=outputs_enc outputs_dec

SRCS:=$(wildcard *.c)
SRCS:=$(filter-out $(TARGET1).c, $(SRCS))
SRCS:=$(filter-out $(TARGET2).c, $(SRCS))
OBJS=$(patsubst %.c,$(OBJDIR)/%.o,$(SRCS))

all: $(TARGET1) $(TARGET2)

debug: CFLAGS += -g -Wall
debug: $(TARGET1) $(TARGET2)

$(OBJDIR)/%.o: %.c | $(OBJDIR) $(OUTDIR)
	$(CC) $(CFLAGS) $? -o $@

$(TARGET1): $(TARGET1).o $(OBJS)
	$(CC) $^ $(LFLAGS) -o $@

$(TARGET2): $(TARGET2).o $(OBJS)
	$(CC) $^ $(LFLAGS) -o $@

$(OBJDIR):
	mkdir -p $@

$(OUTDIR):
	mkdir -p $@

clean:
	rm -r $(OBJDIR)
	rm -r $(OUTDIR)
	rm $(TARGET1)
	rm $(TARGET2)
