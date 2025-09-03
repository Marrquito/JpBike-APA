CXX = g++
CXXFLAGS = -Wall -std=c++17
TARGET = jpBike

# compilar tudo
all: $(TARGET)

$(TARGET): main.cpp
	$(CXX) $(CXXFLAGS) main.cpp -o $(TARGET)

# executar
run: $(TARGET)
	./$(TARGET)

# limpar
clean:
	rm -f $(TARGET)

.PHONY: all run clean
