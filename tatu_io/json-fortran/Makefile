# Build directory path
build = ../../build
# Fortran Compiler
fc = gfortran
# Fortran Compiler flags
flags = -J$(build) -fbacktrace -fno-omit-frame-pointer -O3 -DNDEBUG -O3
# If not exist, create build directory
$(shell mkdir -p $(build))

# json-fortran variable from file
include ./dependencies.make
# Adds build path and .o extension to each one of json-fortran_dependencies
json-fortran_dependencies.o = $(patsubst %, $(build)/%.o, $(json-fortran_dependencies))

# Build the main module object
$(build)/json_module.o: json_module.F90 $(json-fortran_dependencies.o)
	$(fc) $(flags) -c $(<) -o $(@)

# Build dependencies objects
$(build)/%.o: %.F90
	$(fc) $(flags) -c $(<) -o $(@)

clean:
	rm -rf $(build)
