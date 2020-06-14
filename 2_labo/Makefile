CC = g++
CFLAGS = -Wall -O3 --std=c++11

main1.exe : main1.o random.o posizione.o
	$(CC) random.o main1.o posizione.o -o main1.exe
main1.o : main1.cpp
	$(CC) -c main1.cpp -o main1.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
main2.exe : main2.o random.o posizione.o
	$(CC) random.o main2.o posizione.o -o main2.exe
main2.o : main2.cpp
	$(CC) -c main2.cpp -o main2.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
posizione.o : posizione.cpp posizione.h
	$(CC) -c posizione.cpp -o posizione.o $(CFLAGS)
clean1 :
	rm *.o seed.out