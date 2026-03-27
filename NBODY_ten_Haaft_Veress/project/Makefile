# Compiler
CXX = g++

# Compiler flags
CXXFLAGS = -std=c++17 -O2 -Wall -Wextra -Iinclude

# Executable name
TARGET = nbody

# Directories
SRC_DIR = src
INC_DIR = include
BUILD_DIR = build

# Source files
SRCS = $(SRC_DIR)/main.cpp \
       $(SRC_DIR)/Body.cpp \
       $(SRC_DIR)/Euler.cpp \
       $(SRC_DIR)/EulerChromer.cpp \
       $(SRC_DIR)/Hermit.cpp \
       $(SRC_DIR)/Heun.cpp \
       $(SRC_DIR)/IteratedHermit.cpp \
       $(SRC_DIR)/NbodyIntegrator.cpp \
       $(SRC_DIR)/RK4.cpp \
       $(SRC_DIR)/vector.cpp \
       $(SRC_DIR)/VelocityVerlet.cpp

# Object files in build/
OBJS = $(SRCS:$(SRC_DIR)/%.cpp=$(BUILD_DIR)/%.o)

# Default target
all: $(TARGET)

# Link executable
$(TARGET): $(OBJS)
	$(CXX) $(OBJS) -o $(TARGET)

# Compile .cpp into build/.o
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Create build directory if it does not exist
$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

# Clean
clean:
	rm -rf $(BUILD_DIR) $(TARGET)

# Rebuild
rebuild: clean all