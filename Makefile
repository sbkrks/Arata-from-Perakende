CPP_FILES := $(wildcard src/*.cpp)
ExeSuf=9.0
MAINEVENTS    = main.cpp 
MAINEVENTO    = arata$(ExeSuf) 

all:
	g++ -std=c++11 $(CPP_FILES) $(MAINEVENTS)  -o $(MAINEVENTO)
clean:
	rm arata *.txt
