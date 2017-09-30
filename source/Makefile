objects = prop_migra1.o main_migra1.o grid_migra1.o gold_migra1.o funs_migra1.o sim_migra1.o rand.o
macro = -DUNDEBUG
mode = -g -Wall -ffloat-store
compiler = gcc
MACH = `hostname|sed 's/^\([^.]*\).*/\1/'`



migra1	 : $(objects)
	$(compiler) $(mode) -o migra1 $(objects) -lm
	if [ -d ../build ]; then \
	   mv migra1 ../build/migra1-$(MACH) ; \
	else \
	     mv migra1 ../migra1-$(MACH) ;\
	fi
	rm -f $(objects)

prop_migra1.o : prop_migra1.c migra1.h 
	$(compiler) $(mode) $(macro) -c prop_migra1.c 
main_migra1.o : main_migra1.c migra1.h 
	$(compiler) $(mode) $(macro) -c main_migra1.c 
grid_migra1.o : grid_migra1.c migra1.h 
	$(compiler) $(mode) $(macro) -c grid_migra1.c 
gold_migra1.o : gold_migra1.c migra1.h 
	$(compiler) $(mode) $(macro) -c gold_migra1.c 
funs_migra1.o : funs_migra1.c migra1.h 
	$(compiler) $(mode) $(macro) -c funs_migra1.c 
sim_migra1.o : sim_migra1.c migra1.h 
	$(compiler) $(mode) $(macro) -c sim_migra1.c 
rand.o : rand.c
	$(compiler) $(mode) $(macro) -c rand.c 


.PHONY : clean
clean :
	-rm $(objects)
