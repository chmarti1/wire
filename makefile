pre:
	apt install libblas-dev liblapacke-dev

wsolve: wsolve.c wsolve.h
	gcc wsolve.c -lm -llapacke -lpthread -o wsolve
