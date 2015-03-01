all: simulation

simulation: Tensor.o FonctionsAnnexes.o
	g++ -o WL_bin main.cpp Tensor.o FonctionsAnnexes.o -std=c++11 -W -Wall -pedantic
	g++ -o WL_computation_bin main_computation.cpp -std=c++11 -W -Wall -pedantic

Tensor.o: Tensor.cpp
	g++ -o Tensor.o -c Tensor.cpp -std=c++11 -W -Wall -pedantic

FonctionsAnnexes.o: fonctionsAnnexes.cpp
	g++ -o FonctionsAnnexes.o -c fonctionsAnnexes.cpp -std=c++11 -W -Wall -pedantic

clean:
	rm -rf *.o
#	mv Cv* data/
#	mv lng* data/

mkproper: clean
	rm -rf WL_bin WL_computation_bin
