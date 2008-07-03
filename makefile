CC = cc
CFLAGS = -O3 -mpentium

all: testmesh showmesh cube

testmesh: massmesh.o testmesh.o
	$(CC) $(CFLAGS) -o testmesh massmesh.o testmesh.o -lm

showmesh: massmesh.o graphics.o showmesh.o
	$(CC) $(CFLAGS) -o showmesh massmesh.o graphics.o showmesh.o \
	-lglut -lGLU -lGL -lm

cube: massmesh.o graphics.o cube.o
	$(CC) $(CFLAGS) -o cube massmesh.o cube.o graphics.o -lglut -lGLU -lGL -lm

clean:
	-rm *.o testmesh showmesh cube
