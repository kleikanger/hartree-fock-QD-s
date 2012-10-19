L = -std=c++11 -O3

run: hartreeFock.out
	./hartreeFock.out

hartreeFock.out: hartreeFock.o hamiltonianElem.o newmatrix.o libGRIE.o lib.o
	g++ $(L) lib.o newmatrix.o hartreeFock.o libGRIE.o classQDotInteraction.o modGaussTools.o classRadialPotential.o -o hartreeFock.out -llapack

#newmatrix.o: newmatrix.cpp newmatrix.h
#	g++ $(L) -c newmatrix.cpp

hartreeFock.o: hartreeFock.cpp dsyev.h newmatrix.h hamiltonianElem.o libGRIE.o
	g++ $(L) -c hartreeFock.cpp

hamiltonianElem.o: hamiltonianElem.cpp hamiltonianElem.h libGRIE.o
	g++ $(L) -c hamiltonianElem.h hamiltonianElem.cpp

lib.o: lib.cpp lib.h
	g++ $(L) -c lib.cpp

libGRIE.o: libGRIE.cpp libGRIE.h classQDotInteraction.o
	$(CC)  $(L)  -c libGRIE.cpp

classQDotInteraction.o: OpenFCI/classQDotInteraction.cpp OpenFCI/classQDotInteraction.hpp modGaussTools.o classRadialPotential.o 
	$(CC) -c OpenFCI/classQDotInteraction.cpp 

modGaussTools.o: OpenFCI/modGaussTools.cpp OpenFCI/modGaussTools.hpp 
	$(CC) -c OpenFCI/modGaussTools.cpp 

classRadialPotential.o: OpenFCI/classRadialPotential.cpp OpenFCI/classRadialPotential.hpp 
	$(CC) -c OpenFCI/classRadialPotential.cpp 

clear:
	rm *o *out
