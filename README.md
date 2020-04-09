# Kujo
Python code for exciton dynamics in organic single crystals. In development by Igor Koskin (https://github.com/TovCat).

### Current progress
The code still lacks most of the required key features and still is in the early stages of development. Therefore, the validity of the results obtained while using Kujo are not guaranteed. If you happen to stumle upon this code (which you probably shouldn't, but anyway) - use it at your own risk.

### Executing code
Kujo works only in Python 3 environment - currently code is not Python 2 compartible and probably wouldn't in observable future. In order to execute Kujo just run in your prompt:
* python3 main.py <input_file>

whereas <input_file> stands for instruction file, telling the code what to calculate.

Aside from build-in Python libraries, the code uses NumPy and mathplotlib.

### Input file
Syntax of the input file is pretty much straightforward. Instructions are executed line by line:
* (function): (var1) = ...; (var2) = ...; (var3) = ...; (...)

whereas function tells code what to calculate, colon separates function from variables feeded into the function. Variables must contain variable name, equals sign (=) followed by the value of a variable and semicolon to separate one variable from another. Input file is case sensitive. 

### Available functions
#### read_cube
Reads .cube file for further operations. Currently the code is only capable to store one .cube file at a time, so reading another .cube file will erase the previous one.

_Variables_:
* file - name of the .cube file to read. Currently only files in the code directory could be read.

#### read_cif
Reads .cube file for further operations. Currently the code is only capable to store one .cif file at a time, so reading another .cif file will erase the previous one.

_Variables_:
* file - name of the .cif file to read. Currently only files in the code directory could be read.

#### build_cluster
Builds molecular cluster via replication of the primitive cell of preiously read .cif file. Obviously, calling this function before reading .cif file would crush the code.

_Variables_:
* a: number of times to replicate primitive cell in the direction of a-axis. 0 stands for no replications at all.
* b: -//- b-axis.
* c: -//- c-axis.

#### integrate_translated_cube
Calculates exciton coupling between two translated cubes of Transition Density Matrix.

_Variables_:
* vector - coordinates of translation vector in [Angstrom].
* vector-cif - translation vector from .cif file (obviously, available only if you read .cife file before). Possible values are: "a", "b", "c" standing for a, b, c translation vectors respectively.
* multiplier - number (floating or integer) to multiply the vector.

#### set_hard_cutoff
Set cutoff distance after which code would set every coupling equal to 0 no matter what.

_Variables_:
Value - cutoff distance in [Angstrom]

#### set_int_cutoff
Set cutoff distance after which code would switch comptutational model from TDM integration to Extanded Dipole Model in order to save time.

_Variables_:
Value - cutoff distance in [Angstrom]

Currently, not all functions that coded in Kujo are accessible via input file. You may try to use them by importing Kujo libraries directly, however if you are considering to do that, I would urge you to reconsider, since no guaranties are given that they would work the way you want them to work.
