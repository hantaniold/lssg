CC=mpicc
# GA libraries are installed in /usr/local/lib
# GA headers are installed in /usr/local/include
CPPFLAGS=-I/usr/local/include
LDFLAGS=-L/usr/local/lib
CFLAGS=-O2 -funroll-loops
LIBS=-lga -L/usr/lib/gcc/x86_64-linux-gnu/4.6.1 -L/usr/lib/gcc/x86_64-linux-gnu/4.6.1/../../../x86_64-linux-gnu -L/usr/lib/gcc/x86_64-linux-gnu/4.6.1/../../../../lib -L/lib/x86_64-linux-gnu -L/lib/../lib -L/usr/lib/x86_64-linux-gnu -L/usr/lib/../lib -L/usr/lib/gcc/x86_64-linux-gnu/4.6.1/../../.. -lgfortran -lm -lquadmath
FC = mpif90
CPPFLAGS += -DUSE_MPI -DMPI
LINK = $(FC)
LOADER_OPTS = -g
.SUFFIXES: .c .o .h .x
.PHONY: all
all: matrix-c.x  ga-delta.x
.c.o:
	$(CC) $(CFLAGS) $(CPPFLAGS) -c -o $@ $<
.o.x:
	$(LINK) $(LOADER_OPTS) $(LDFLAGS) -o $@ $< $(LIBS)
clean:
	$(RM) *.o *.x
