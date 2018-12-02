# Build directory path, where to put .o and .mod
build = ./build

# json-fortran dependencies
json-fortran = json_module json_kinds json_parameters json_value_module json_string_utilities json_file_module
# json_io dependencies
json_io = cli file input_types json_io $(json-fortran)
# Adds json_io build path and .o extension to each name of json_io variable
json_io.o = $(patsubst %, $(build)/%.o, $(json_io))

# If not exist, create build directory
$(shell mkdir -p $(build))
# Electromagnetics dipoles dependencies
dependencies = Anderson computationalstuff filtros escolhadofiltro dehx dehy dev dmhx dmhy dmv csv_file main1D
# Adds build path and .o extension to each one of json_io_files
dependencies.o = $(patsubst %, $(build)/%.o, $(dependencies))


# Target to create executable binary
mod1d.x: $(json_io.o) $(dependencies.o)
# $(@) represents the current target, in this case: $(@) = mod1d.x
# $(^) represents all dependencies of the current target, in this case: all .o files
# -J specifies where to search for .mod files for compiled modules
	gfortran -J$(build) -o $(@) $(^)

# Runs the first target of Makefile in ./json_io/json_io directory
# Compile json_io module and all its dependencies
$(build)/%.o: json_io/source/%.f08
	$(MAKE) -C json_io

$(build)/Anderson.o: Anderson.for
	gfortran -J$(build) -std=legacy -c $(<) -o $(@)

$(build)/%.o: %.f08
# $(<) represents the first dependency of the current target, in this case: $(<) = %.f08
	gfortran -J$(build) -std=f2008 -pedantic -c $(<) -o $(@)

clean:
	rm -rf $(build)
	$(MAKE) -C json_io clean
