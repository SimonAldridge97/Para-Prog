CXX = g++
CXXFLAGS = -O2 -std=c++17

TARGET = mergesort
SRC = mergesort.cpp

all: $(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(SRC)

clean:
	rm -f $(TARGET) results.txt benchmark_plot.png
