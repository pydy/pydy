all : rattleback_sim rattleback_stability

rattleback.o : rattleback.cc
	g++ -Wall -O3 -c rattleback.cc -I/home/luke/usr/include -I/home/luke/usr/include/eigen3

rattleback_sim.o : rattleback_sim.cc
	g++ -Wall -O3 -c rattleback_sim.cc -I/home/luke/usr/include

rattleback_stability.o : rattleback_stability.cc
	g++ -Wall -O3 -c rattleback_stability.cc -I/home/luke/usr/include -I/home/luke/usr/include/eigen3

rattleback_sim : rattleback_sim.o rattleback.o
	g++ -Wall -O3 -o rattleback_sim rattleback_sim.o rattleback.o -lgsl -lgslcblas

rattleback_stability : rattleback_stability.o rattleback.o
	g++ -Wall -O3 -o rattleback_stability rattleback_stability.o rattleback.o -lgsl -lgslcblas

simulation.dat : rattleback_sim
	./rattleback_sim

reproduceKane1982.o : reproduceKane1982.cc
	g++ -Wall -O3 -c reproduceKane1982.cc -I/home/luke/usr/include

reproduceKane1982_sim.o : reproduceKane1982_sim.cc
	g++ -Wall -O3 -c reproduceKane1982_sim.cc -I/home/luke/usr/include

reproduceKane1982_sim : reproduceKane1982_sim.o reproduceKane1982.o
	g++ -Wall -O3 -o reproduceKane1982_sim reproduceKane1982_sim.o reproduceKane1982.o -lgsl -lgslcblas

clean :
	rm -rf *.o rattleback_sim rattleback_stability reproduceKane1982_sim simulation.dat

.PHONY : all clean
