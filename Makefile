# Build directory path, where to put .o and .mod
build = ./build
# Fortran Compiler
fc = gfortran -J$(build) -std=f2008
# Flags to Fortran Compiler
# -J specifies where to put .mod files for compiled modules
development_flags = -pedantic -Wall -Wextra -Wimplicit-interface -fPIC -fmax-errors=1 -g -fcheck=all -fbacktrace
production_flags  = -pedantic -Wall -Wextra -Wimplicit-interface -fPIC -Werror \
 -fmax-errors=1 -O3 -march=native -ffast-math -funroll-loops -static-libgfortran

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
include tatu_io/json-fortran/dependencies.make
# tatu_io dependencies from file
include tatu_io/dependencies.make
# Adds build path and .o extension to each one of json-fortran dependencies
json-fortran.o = $(patsubst %, $(build)/%.o, $(json-fortran))
# Adds build path and .o extension to each one of json_io dependencies
tatu_io.o = $(patsubst %, $(build)/%.o, $(tatu_io))
#--------------------------------------------------------------------

# TATU
#--------------------------------------------------------------------
# tatu dependencies from file
include ./dependencies.make
# Adds build path and .o extension to each one of tatu dependencies
tatu.o = $(patsubst %, $(build)/%.o, $(tatu))
#--------------------------------------------------------------------

# GNU Make Target-specific Variable Values
development: flags = $(development_flags)
production: flags = $(production_flags)

# Development target (default)
development: tatu

# Production target
production: tatu

# Target to create executable binary
tatu: $(tatu_io.o) $(json-fortran.o) $(clifor.o) $(tatu.o)
# $(@) represents the current target, in this case: $(@) = tatu
# $(^) represents all dependencies of the current target, in this case: all .o files
# -J specifies where to search for .mod files for compiled modules
	$(fc) $(flags) -o $(@) $(^)

# Compile clifor module and all its dependencies
$(build)/%.o: clifor/source/%.f08
	$(MAKE) -C clifor

# Runs the first target of Makefile in ./tatu_io/tatu_io directory
# Compile tatu_io module and all its dependencies (json-fortran included)
$(build)/%.o: tatu_io/source/%.f08
	$(MAKE) -C tatu_io

$(build)/Anderson.o: Anderson.for
	$(fc) -std=legacy -c $(<) -o $(@)

$(build)/%.o: %.f08
# $(<) represents the first dependency of the current target, in this case: $(<) = %.f08
	$(fc) $(flags) -c $(<) -o $(@)

# Compile files inside models subdirectory
$(build)/%.o: models/%.f08
	$(fc) $(flags) -c $(<) -o $(@)

clean:
	rm -rf $(build)
