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
flags = -J$(build) -std=f2008 -pedantic -Wall -Wextra -Wimplicit-interface -fPIC -fmax-errors=1 -g -fcheck=all -fbacktrace

# CLIFor
#--------------------------------------------------------------------
# clifor dependencies from file
include ../clifor/dependencies.make
# Adds build path and .o extension to each one of clifor dependencies
clifor.o = $(patsubst %, $(build)/%.o, $(clifor))
#--------------------------------------------------------------------

# JSON_IO
#--------------------------------------------------------------------
# json-fortran dependencies from file
include ../tatu_io/json-fortran/dependencies.make
# tatu_io dependencies from file
include ../tatu_io/dependencies.make
# Adds build path and .o extension to each one of json-fortran dependencies
json-fortran.o = $(patsubst %, $(build)/%.o, $(json-fortran))
# Adds build path and .o extension to each one of json_io dependencies
tatu_io.o = $(patsubst %, $(build)/%.o, $(tatu_io))
#--------------------------------------------------------------------

# Get all tatu dependencies in a $(tatu) variable from file
include ../dependencies.make
# Object files relative to each source file in $(tatu)
tatu.src = $(filter-out tatu, $(tatu))
tatu.obj = $(patsubst %, $(build)/%.o, $(tatu.src))

# All files matching test_*.$(ext) pattern in test directory
test.src = $(wildcard test_*.$(ext))
# Object files relative to each test file in $(test_src)
test.obj = $(patsubst test_%.$(ext), $(build)/test_%.o, $(test.src))

# If not exist, create build directory
$(shell mkdir -p $(build))


# Run tests
# 1>data/output.txt 2>data/error.txt
default: $(executable)
	./$(executable)
# -i ./data/input.json -o ./data/output

# Create executable binary
$(executable): $(json-fortran.o) $(tatu_io.o) $(clifor.o) $(tatu.obj) $(build)/fruit.o $(test.obj) $(build)/$(driver).o
	$(fc) $(flags) -o $(@) $(^)

# Compile json-fortran files
$(build)/%.o: ../tatu_io/json-fortran/%.F90
	$(fc) -J$(build) -fbacktrace -fno-omit-frame-pointer -O3 -DNDEBUG -O3 -c $(<) -o $(@)

# Compile tatu_io files
$(build)/%.o: ../tatu_io/source/%.f08
	$(fc) $(flags) -c $(<) -o $(@)

# Compile clifor files
$(build)/%.o: ../clifor/source/%.f08
	$(fc) $(flags) -c $(<) -o $(@)

# Compile models files
$(build)/%.o: $(src)/models/%.$(ext)
	$(fc) $(flags) -c $(<) -o $(@)

# Compile tatu files
$(build)/%.o: $(src)/%.$(ext)
	$(fc) $(flags) -c $(<) -o $(@)

$(build)/Anderson.o: $(src)/Anderson.for
	$(fc) -std=legacy -c $(<) -o $(@)

# Compile fruit test framework
$(build)/fruit.o: ./fruit/fruit.f90
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
