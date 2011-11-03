all : rattleback_sim

rattleback.o : rattleback.cc
	g++ -Wall -O3 -c rattleback.cc -I/home/luke/usr/include -I/home/luke/usr/include/eigen3

rattleback_sim.o : rattleback_sim.cc
	g++ -Wall -O3 -c rattleback_sim.cc -I/home/luke/usr/include

rattleback_sim : rattleback_sim.o rattleback.o
	g++ -Wall -O3 -o rattleback_sim rattleback_sim.o rattleback.o -lgsl -lgslcblas

clean :
	rm -rf *.o rattleback_sim datafile.dat

.PHONY : all clean
