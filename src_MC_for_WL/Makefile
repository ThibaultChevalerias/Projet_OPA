all: simulation

simulation: Tensor.o FonctionsAnnexes.o
	g++ -o classic_MC_bin main.cpp Tensor.o FonctionsAnnexes.o -std=c++11 -W -Wall -pedantic

Tensor.o: Tensor.cpp
	g++ -o Tensor.o -c Tensor.cpp -std=c++11 -W -Wall -pedantic

FonctionsAnnexes.o: fonctionsAnnexes.cpp
	g++ -o FonctionsAnnexes.o -c fonctionsAnnexes.cpp -std=c++11 -W -Wall -pedantic

clean:
	rm -rf *.o

mkproper: clean
	rm -rf classic_MC_bin
