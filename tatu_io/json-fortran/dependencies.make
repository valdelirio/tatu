# Makefile to be included in another ones to variables sharing

# json-fortran dependencies
json-fortran_dependencies = json_kinds json_parameters json_string_utilities json_value_module json_file_module
# All files of json-fortran
json-fortran = $(json-fortran_dependencies) json_module
