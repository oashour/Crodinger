CC=gcc
CFLAGS=-Wall -std=c99
INCLUDES=-Ilib
LDFLAGS=-lfftw3 -lm
DIR = lib
SRCS = nlse.c $(DIR)/lib.c
OBJS = $(SRCS:.c=.o)
DEPS = lib.h
MAIN = nlse

.PHONY: clean

all: $(MAIN)
	@echo NLSE solver has been compiled

$(MAIN): $(OBJS) 
	$(CC) -o $@ $^ $(INCLUDES) $(CFLAGS) $(LDFLAGS) 

.c.o: $(SRCS) $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS) $(LDFLAGS) $(INCLUDES)

clean:
	$(RM) *.o $(DIR)/*.o *~ $(MAIN)
