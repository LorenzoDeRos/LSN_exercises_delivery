CC = g++
CFLAGS = -Wall -O3 --std=c++11

tsp_metropolis.exe : tsp_metropolis.o random.o
	$(CC) random.o tsp_metropolis.o -o tsp_metropolis.exe
tsp_metropolis.o : tsp_metropolis.cpp genetic.h
	$(CC) -c tsp_metropolis.cpp -o tsp_metropolis.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
tsp_parallel_genetic.exe : tsp_parallel_genetic.o random.o
	mpicxx random.o tsp_parallel_genetic.o -o tsp_parallel_genetic.exe
tsp_parallel_genetic.o : tsp_parallel_genetic.cpp genetic.h
	mpicxx -c tsp_parallel_genetic.cpp -o tsp_parallel_genetic.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o *.out
