CC = g++
CFLAGS = -Wall -O3 --std=c++11
LIBS:=`root-config --libs`
INCS:=`root-config --cflags`

main.exe : main.o random.o posizione.o
	$(CC) random.o main.o posizione.o -o main.exe
main.o : main.cpp
	$(CC) -c main.cpp -o main.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
posizione.o : posizione.cpp posizione.h
	$(CC) -c posizione.cpp -o posizione.o $(CFLAGS)
%.o : %.cpp
	g++ -Wall -c $< ${INCS}
qmc1d: qmc1d.o
	g++ -O3 -Wall -o $@ $^ ${LIBS}
clean:
	rm *.o seed.out potential.dat kinetic.dat probability.dat
