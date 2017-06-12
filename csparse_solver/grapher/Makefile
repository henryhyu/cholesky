CC=gcc
CFLAGS=-I.
DEPS=Solver.h
OBJ=Solver.o csparse.o mmio.o 
LIBS=-lm

CFLAGS += $(shell pkg-config --cflags json-c)
LIBS += $(shell pkg-config --libs json-c)

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS) $(LIBS)

solve: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

clean:
	rm -f *.o solve