OBJECTS = conj_grad.o soft_sphere.o main.o 
CXXFLAGS =  -O3 
min.out : $(OBJECTS) 
	g++ -o min.out $(OBJECTS)
conj_grad.o : conj_grad.cpp conj_grad.h
soft_sphere.o : soft_sphere.cpp soft_sphere.h


