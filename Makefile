CC = g++

CXXFLAGS = -g -Wall -I/usr/include/eigen3 -I/usr/include/mpi #-D_CC_OVERLAP 

LDFLAGS = -lm 

GFOBJS = parameters.o mb_self_energy.o interacting_gf.o current.o leads_self_energy.o sigma_2.o##impurity_solver.o green_function.o

EXECS = multi_orbital_sigma_2

##i think the problem is that i dont have explicitly the dmft depends on parameters.

all: $(EXECS)

multi_orbital_sigma_2: main.o  $(GFOBJS) 
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

tags:
	etags *.c *.h

.PHONY: clean tags tests

clean:
	$(RM) *.o $(EXECS) $(TESTS) TAGS tags
