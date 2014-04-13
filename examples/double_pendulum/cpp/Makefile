all : double_pendulum_sim

double_pendulum_ode.o : double_pendulum.cc
	g++ -Wall -c double_pendulum.cc

double_pendulum_sim : double_pendulum_sim.cc double_pendulum.o
	g++ -Wall double_pendulum_sim.cc double_pendulum.o -lgsl -lgslcblas -o double_pendulum_sim

clean :
	rm -rf *.o double_pendulum_sim

.PHONY : all clean
