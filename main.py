import excitons
import reader_cif
import reader_orca
import reader_cube
import numpy as np
import linker
from os import getcwd
import kujo_io


instructions, options = kujo_io.read_input("input.txt")