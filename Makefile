# Build directory path, where to put .o and .mod
build = ./build
# Fortran Compiler
fc = gfortran
# Flags to Fortran Compiler
# -J specifies where to put .mod files for compiled modules
development_flags = -J$(build) -std=f2008 -pedantic -Wall -Wextra -Wimplicit-interface -fPIC -fmax-errors=1 -g -fcheck=all -fbacktrace
production_flags  = -J$(build) -std=f2008 -pedantic -Wall -Wextra -Wimplicit-interface -fPIC -Werror \
 -fmax-errors=1 -O3 -march=native -ffast-math -funroll-loops -static-libgfortran
flags = $(production_flags)

# If not exist, create build directory
$(shell mkdir -p $(build))

# CLIFor
#--------------------------------------------------------------------
# clifor dependencies from file
include clifor/dependencies.make
# Adds build path and .o extension to each one of clifor dependencies
clifor.o = $(patsubst %, $(build)/%.o, $(clifor))
#--------------------------------------------------------------------

# JSON_IO
#--------------------------------------------------------------------
# json-fortran dependencies from file
include json_io/json-fortran/dependencies.make
# json_io dependencies from file
include json_io/dependencies.make
# Adds build path and .o extension to each one of json-fortran dependencies
json-fortran.o = $(patsubst %, $(build)/%.o, $(json-fortran))
# Adds build path and .o extension to each one of json_io dependencies
json_io.o = $(patsubst %, $(build)/%.o, $(json_io))
#--------------------------------------------------------------------

# ELECTROMAGNETICS DIPOLES
#--------------------------------------------------------------------
# dipoles dependencies from file
include ./dependencies.make
# Adds build path and .o extension to each one of dipoles dependencies
dipoles.o = $(patsubst %, $(build)/%.o, $(dipoles))
#--------------------------------------------------------------------


# Target to create executable binary
tatu: $(json_io.o) $(json-fortran.o) $(clifor.o) $(dipoles.o)
# $(@) represents the current target, in this case: $(@) = main1D.x
# $(^) represents all dependencies of the current target, in this case: all .o files
# -J specifies where to search for .mod files for compiled modules
	$(fc) $(flags) -o $(@) $(^)

# Compile clifor module and all its dependencies
$(build)/%.o: clifor/source/%.f08
	$(MAKE) -C clifor

# Runs the first target of Makefile in ./json_io/json_io directory
# Compile json_io module and all its dependencies (json-fortran included)
$(build)/%.o: json_io/source/%.f08
	$(MAKE) -C json_io

$(build)/Anderson.o: Anderson.for
	$(fc) -J$(build) -std=legacy -c $(<) -o $(@)

$(build)/%.o: %.f08
# $(<) represents the first dependency of the current target, in this case: $(<) = %.f08
	$(fc) $(flags) -c $(<) -o $(@)

# Compile files inside models subdirectory
$(build)/%.o: models/%.f08
	$(fc) $(flags) -c $(<) -o $(@)

clean:
	rm -rf $(build)
