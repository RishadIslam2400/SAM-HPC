# Compiler
NVCC = nvcc

# Compiler flags
CFLAGS = -O2 -std=c++11

# Executable name
TARGET = main3

# Source files
SRCS = main3.cu

# Object files
OBJS = $(SRCS:.cu=.o)

# Default rule
all: $(TARGET)

# Rule to link the target executable
$(TARGET): $(OBJS)
	$(NVCC) $(CFLAGS) -o $(TARGET) $(OBJS)

# Rule to compile source files into object files
%.o: %.cu
	$(NVCC) $(CFLAGS) -c $< -o $@

# Clean rule
clean:
	rm -f $(OBJS) $(TARGET)

# Phony targets
.PHONY: all clean