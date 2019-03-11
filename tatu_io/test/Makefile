# Source extension
ext = f08
# Test driver basename
driver = driver
# Test binary executable
executable = driver.bin
# Source path
src = ../source
# Build path
build = ./build
# Fortran Compiler
fc = gfortran
# Fortran Compiler flags
flags = -J$(build) -std=f2008 -pedantic

# Get all json-fortran dependencies in a $(json-fortran) variable from file
include ../json-fortran/dependencies.make
# Adds build path and .o extension to each one of json-fortran dependencies
json-fortran.o = $(patsubst %, $(build)/%.o, $(json-fortran))

# Get all tatu_io dependencies in a $(tatu_io) variable from file
include ../dependencies.make
# Adds build path and .o extension to each one of tatu_io dependencies
tatu_io.o = $(patsubst %, $(build)/%.o, $(tatu_io))

# All files matching test_*.$(ext) pattern in test directory
test_src = $(wildcard test_*.$(ext))
# Object files relative to each test file in $(test_src)
test_obj = $(patsubst test_%.$(ext), $(build)/test_%.o, $(test_src))

# If not exist, create build directory
$(shell mkdir -p $(build))


# Run tests
default: $(executable)
	./$(executable) -i ./data/input.json -o ./data/output.json

# Create executable binary
$(executable): $(json-fortran.o) $(tatu_io.o) $(build)/fruit.o $(test_obj) $(build)/$(driver).o
	$(fc) $(flags) -o $(@) $(^)

# Compile json-fortran files
$(build)/%.o: ../json-fortran/%.F90
	$(fc) $(flags) -c $(<) -o $(@)

# Compile tatu_io files
$(build)/%.o: $(src)/%.$(ext)
	$(fc) $(flags) -c $(<) -o $(@)

# Compile fruit test framework
$(build)/fruit.o: ../fruit/fruit.f90
	$(fc) $(flags) -c $(<) -o $(@)

# Compile test files
$(build)/test_%.o: test_%.$(ext)
	$(fc) $(flags) -c $(<) -o $(@)

# Compile test driver
$(build)/$(driver).o: $(driver).$(ext)
	$(fc) $(flags) -c $(<) -o $(@)

.PHONY: clean cleanall

clean:
	rm -rf $(build)
