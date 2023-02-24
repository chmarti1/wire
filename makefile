pre:
	apt install libblas-dev liblapacke-dev

wireslice: wireslice.c wireslice.h
	gcc wireslice.c -lm -llapacke -lpthread -o wireslice
