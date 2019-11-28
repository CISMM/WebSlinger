CC = g++
CFLAGS = -O3 

all: meshmaker testmesh showmesh cube webslinger

meshmaker: meshmaker.o
	$(CC) $(CFLAGS) -o meshmaker meshmaker.o -lm

testmesh: massmesh.o testmesh.o
	$(CC) $(CFLAGS) -o testmesh massmesh.o testmesh.o -lm

showmesh: massmesh.o graphics.o showmesh.o
	$(CC) $(CFLAGS) -o showmesh massmesh.o graphics.o showmesh.o \
	-lglut -lGLU -lGL -lm

cube: massmesh.o graphics.o cube.o
	$(CC) $(CFLAGS) -o cube massmesh.o cube.o graphics.o -lglut -lGLU -lGL -lm

webslinger: graphics.o massmesh.o webslinger.o
	$(CC) $(CFLAGS) -o webslinger graphics.o massmesh.o webslinger.o -lglut -lGLU -lGL -lm

clean:
	-rm *.o testmesh showmesh cube webslinger
