all: linsker

linsker: Main.o Network.o
	g++ Main.o Network.o -o linsker

Main.o: Main.cpp
	g++ -c Main.cpp

Network.o: Network.cpp
	g++ -c Network.cpp

clean:
	rm -rf *.o linsker