all : double_pendulum_sim

double_pendulum_ode.o : double_pendulum_ode.cc
	g++ -Wall -c double_pendulum_ode.cc

double_pendulum_sim : double_pendulum_sim.cc double_pendulum_ode.o
	g++ -Wall double_pendulum_sim.cc double_pendulum_ode.o  -o double_pendulum_sim -lgsl -lgslcblas

clean :
	rm -rf double_pendulum_ode.o double_pendulum_sim

.PHONY : all clean
