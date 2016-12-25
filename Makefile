CPP_FILES := $(wildcard src/*.cpp)
MAINEVENTS    = main.cpp 
MAINEVENTO         = arata$(ExeSuf) 

all:
	g++ -std=c++11 $(CPP_FILES) $(MAINEVENTS)  -o $(MAINEVENTO)
compile:
