import os
import tkinter as tk 
from tkinter.filedialog import askopenfilename,  askdirectory
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Polypeptide import PPBuilder
import ntpath
import warnings
from Bio import BiopythonWarning
from get_raw_distances import *


model = "/home/rotan/Desktop/Objective2/s_test.pdb"
parser = PDBParser(PERMISSIVE=1)
structure_id = ntpath.basename(model).split('.')[0]
filename = ntpath.basename(model)
structure = parser.get_structure(structure_id, filename)
ppb=PPBuilder()
ppb.build_peptides(structure)
start_pos = []
for pp in ppb.build_peptides(structure):
    start_pos = int((pp[0].get_id()[1]))
    sequence = pp.get_sequence()
cmd.load(model)	
		
cmd.select("residue829", "resi " + str(829))
for i in range(len(sequence)):
    print(('resiude' + str(i+start_pos)))
    residue_a = cmd.select("residue_a", "resi " + str(start_pos + i))
    cmd.distance(('resiude' + str(i+start_pos)), "residue_a", "residue829", "4", mode="2")
    D = get_raw_distances(('resiude' + str(i+start_pos)))
    print("number of h-bonds:", len(D))
