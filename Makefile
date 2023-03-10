.PHONY: default clean

SOURCE	= $(filter-out collapse.c, $(wildcard *.c))
EXEC	= $(SOURCE:.c=)

CC	= gcc
CFLAGS	= -fopenmp -O3
LDFLAGS	= -fopenmp
LDLIBS	= -lm

default: $(EXEC)

clean:
	rm -f $(EXEC)
