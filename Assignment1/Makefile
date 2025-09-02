# Compiler and flags
CXX = g++
CXXFLAGS = -O2 -std=c++17 -Wall

# Target executable
TARGET = mergesort

# Source files
SRC = mergesort.cpp

# Default rule
all: $(TARGET)

# Build rule
$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(SRC)

# Clean rule
clean:
	rm -f $(TARGET) results.txt benchmark_plot.png

.PHONY: all clean
