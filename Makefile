CXX = g++
CXXFLAGS = -Wall -std=c++17 -Iinclude
TARGET = jpBike

# arquivos fonte
SOURCES = main.cpp src/Instance.cpp
OBJECTS = $(SOURCES:.cpp=.o)

# compilar tudo
all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(OBJECTS) -o $(TARGET)

# compilar arquivos .cpp em .o
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# executar
run: $(TARGET)
	./$(TARGET)

# limpar
clean:
	rm -f $(TARGET) $(OBJECTS)

.PHONY: all run clean
