import ovito
import ovito.pipeline
import scipy
import ase
import math
import scipy
import numpy as np
from ovito.io import import_file, export_file
from ovito.modifiers import FreezePropertyModifier ,CombineDatasetsModifier, ClusterAnalysisModifier, DeleteSelectedModifier, ExpressionSelectionModifier, WignerSeitzAnalysisModifier, InvertSelectionModifier
from ase.visualize import view
from ase.build import fcc110, fcc111, diamond111, diamond100, surface, bulk
from ase import Atoms
from ase.cluster import Octahedron
from ovito.io.ase import ovito_to_ase, ase_to_ovito
from ovito.pipeline import Pipeline, StaticSource
from ase.io import write, read
from ase.collections import g2




lattice_constant = 5.468
size = 5



plane = diamond111('Si',size=(int(size/2),int(size/2), size), periodic=True , a = lattice_constant)
#plane = surface('Si', (1,1,0), 9, periodic= True)
#plane.center(vacuum =, axis=2)
#print(n)
plane = bulk('NaCl', 'fcc').repeat((size, size, size))
#plane.center(vacuum = 10, axis=2)
n = len(plane.numbers[...])
# Visualize the structure
print(n)
view(plane)

"""
can just change this to write function
ase_data = ase_to_ovito(plane)
pipeline = Pipeline(source = StaticSource(data=ase_data))
data = pipeline.compute()

export_file(data, "ase_perfect_cube2", 'xyz',columns =
                ["Particle Identifier", "Particle Type",  "Position.X", "Position.Y", "Position.Z"])
"""

#pipeline = import_file("400eV test files/400eV_test2_cluster_1_extract_sample3")
#data = pipeline.compute()
#ase_atoms = ovito_to_ase(data)


#more_atoms = read("400eV test files/400eV_test2_cluster_1_extract_sample1")
#ase_atoms = more_atoms
#ase_atoms += plane
#view(ase_atoms)