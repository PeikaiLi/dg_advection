# Compiler to use
CXX = g++

# Compiler flags:
# -O2       → optimization level 2
# -std=c++17 → use C++17 standard
CXXFLAGS = -O2 -std=c++17

# Linker flags:
# -lopenblas → link against OpenBLAS library
LDFLAGS = -lopenblas

# Output program name
TARGET = test.exe

# Source file(s)
SRC = lapack_test.cpp

# Default rule: build the executable
$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(SRC) $(LDFLAGS)

# Run the program (compile first if needed)
run: $(TARGET)
	./$(TARGET)

# Clean up generated files
clean:
	rm -f $(TARGET)
