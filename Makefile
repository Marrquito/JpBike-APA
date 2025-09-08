CXX = g++
CXXFLAGS = -Wall -std=c++17 -Iinclude -g -O0 -fsanitize=address -fsanitize=undefined
TARGET = jpBike

# arquivos fonte
SOURCES = main.cpp src/Instance.cpp src/Solver.cpp

# compilar tudo
all: $(TARGET)

$(TARGET): $(SOURCES)
	$(CXX) $(CXXFLAGS) $(SOURCES) -o $(TARGET)

# executar
run: $(TARGET)
	./$(TARGET)

# limpar
clean:
	rm -f $(TARGET)

.PHONY: all run clean
