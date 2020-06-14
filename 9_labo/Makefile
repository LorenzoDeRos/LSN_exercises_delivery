CC = g++
CFLAGS = -Wall -O3 --std=c++11

tsp_genetic.exe : tsp_genetic.o random.o
	$(CC) random.o tsp_genetic.o -o tsp_genetic.exe
tsp_genetic.o : tsp_genetic.cpp genetic.h
	$(CC) -c tsp_genetic.cpp -o tsp_genetic.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o *.out