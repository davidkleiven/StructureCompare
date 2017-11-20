import sys
sys.path.append("/home/dkleiven/Documents/StructureCompare")
sys.path.append("/home/davidkl/Documents/StructureCompare")
from ase.build import bulk
import pystructcomp_cpp as pycmp
import numpy as np

def main():
    atoms = bulk( "Al" )
    atoms = atoms*(4,4,4)

    symbols = [atom.symbol for atom in atoms]
    pos = atoms.get_positions()+6.0
    cell = atoms.get_cell().T
    atoms.set_positions(pos)
    atoms.wrap()

    pycmp.test_atoms( symbols, pos, cell )

    print ("Python output:")
    print ("Cell:")
    print (cell)
    print ("Volume:")
    print (np.linalg.det(cell))
    print ("Inverse cell:")
    print (np.linalg.inv(cell))
    print (atoms.get_positions())

if __name__ == "__main__":
    main()
