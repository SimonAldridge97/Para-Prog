#compiler
CXX = g++
# comp flags
CXXFLAGS = -std=c++17 -Wall -O2

# exe name
TARGET = nbody_sim

#source files
SRCS = main.cpp

#obj files
OBJS = $(SRCS:.cpp=.o)

#default rule
all: $(TARGET)

#build exe
$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJS)

# compile source files into objects
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# clean
clean:
	rm -f $(OBJS) $(TARGET)

#run
run: $(TARGET)
	./$(TARGET)




