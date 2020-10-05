CC = gcc
CXX = g++ -std=c++11
DEBUG = 
CXXFLAGS = $(shell fltk-config --use-gl --use-images --cxxflags ) -I.
LDFLAGS = $(shell fltk-config --use-gl --use-images --ldflags )
LDSTATIC = $(shell fltk-config --use-gl --use-images --ldstaticflags )

TARGET = app
OBJS = main.o libtarga.o imageWidget.o scriptHandler.o targaImage.o
SRCS = main.cpp libtarga.c imageWidget.cpp scriptHandler.cpp targaImage.cpp

.SUFFIXES: .o .cpp
*.o: *.cpp
	$(CXX) $(CXXFLAGS) $(DEBUG) -c $< -O3

all: $(TARGET)
	$(CXX) -o $(TARGET) $(OBJS) $(LDSTATIC)

$(TARGET): $(OBJS)

.PHONY: clean
clean:
	rm -f *.o 2> /dev/null
	rm -f $(TARGET) 2> /dev/null

# run script command
# make clean && make all && ./app