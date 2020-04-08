# Kujo
Python code for exciton dynamics in organic single crystals. In development by Igor Koskin (https://github.com/TovCat).

### Current progress
The code still lacks most of the required key features and still in the early stages of development. Therefore, the validity of the results obtained while using Kujo are not guaranteed. If you happen to stumle upon this code (which you probably shouldn't, but anyway) - use it at your own risk.

### Executing code
Kujo works only in Python 3 environment - currently code is not Python 2 compartible and probably wouldn't in observable future. In order to execute Kujo just run in your prompt:
* python3 main.py <input_file>
whereas <input_file> stands for instruction file, telling the code what to calculate.

### Input file
Syntax of the input file is pretty much straightforward. Instructions are executed line by line:
* (function): (var1) = ...; (var2) = ...; (var3) = ...; (...)

whereas function tells code what to calculate, colon separates function from variables feeded into the function. Variables must contain variable name, equals sign (=) followed by the value of a variable and semicolon to separate one variable from another. Input file is case sensitive. 

### Available functions
#### read_cube
Reads .cube file for further operations. Currently the code is only capable to store one .cube file at a time, so reading another .cube file will erase the previous one.
Variables:
* file = <...> - name of the .cube file to read. Currently only files in the code directory could be read.

#### read_cif
Reads .cube file for further operations. Currently the code is only capable to store one .cif file at a time, so reading another .cif file will erase the previous one.
Variables:
* file = <...> - name of the .cif file to read. Currently only files in the code directory could be read.

#### integrate_translated_cube
Calculates exciton coupling between two translated cubes of Transition Density Matrix.
Variables:
* vector: coordinates of translation vector in [Angstrom].
* vector-cif: translation vector from .cif file (obviously, available only if you read .cife file before). Possible values are: "a", "b", "c" standing for a, b, c translation vectors respectively.
* multiplier: number (floating or integer) to multiply the vector.

Currently, not all functions that coded in Kujo are accessible via input file. You may try to use them by importing Kujo libraries directly, however if considering to do that, I would urge you to reconsider your decision, since not guaranties are given that they would work the way you want them to work.
