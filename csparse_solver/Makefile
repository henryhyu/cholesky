CC=gcc
CFLAGS=-I.
DEPS=Solver.h
OBJ=Solver.o csparse.o mmio.o 
LIBS=-lm

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

solve: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

clean:
	rm -f *.o solve