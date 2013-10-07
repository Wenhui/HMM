hmm: hmm.o main.cpp
	g++ -g -Wall -o hmm main.cpp hmm.o
hmm.o: hmm.h hmm.cpp
	g++ -g -Wall -c hmm.cpp
