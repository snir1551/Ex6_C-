#!make -f

CXX=clang++-9 
CXXFLAGS=-std=c++2a

HEADERS := solver.hpp $(wildcard *.h*)
STUDENT_SOURCES := $(filter-out $(wildcard Test*.cpp), $(wildcard *.cpp))
STUDENT_OBJECTS := $(subst .cpp,.o,$(STUDENT_SOURCES))

run: test
	./$^

test: TestRunner.o Test_iris.o Test_shahar.o $(STUDENT_OBJECTS)
	$(CXX) $(CXXFLAGS) $^ -o test

%.o: %.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) --compile $< -o $@

clean:
	rm -f *.o test
