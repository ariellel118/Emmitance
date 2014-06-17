CC = g++
LDFLAGS = `pkg-config --libs ibsimu-1.0.5new_solver`
CXXFLAGS = -Wall -g `pkg-config --cflags ibsimu-1.0.5new_solver`

vlasovE: vlasovE.o
	$(CC) -o vlasovE vlasovE.o $(LDFLAGS)

vlasovE.o: vlasovE.cpp
	$(CC) -c -o vlasovE.o vlasovE.cpp $(CXXFLAGS)

clean:
	$(RM) *~ *.o vlasovE


