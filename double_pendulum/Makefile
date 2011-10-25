all : double_pendulum_sim

double_pendulum_ode.o : double_pendulum_ode.cc
	g++ -Wall -c double_pendulum_ode.cc

double_pendulum_sim : double_pendulum_sim.cc double_pendulum_ode.o
	g++ -Wall -c double_pendulum_sim.cc
	g++ -Wall double_pendulum_ode.o double_pendulum_sim.o -lgsl -lgslcblas -o double_pendulum_sim

clean :
	rm -rf *.o double_pendulum_sim

.PHONY : all clean
