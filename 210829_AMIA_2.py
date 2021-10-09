# This script is developed for the fufuillment for Masters at the South African National Bioinformatics Institute at the University of the Western Cape
# Project Title: "Automated computational workflow to prioritize potential resistance variants identified in HIV Integrase Subtype C and CRF02_AG" 
# Currently the project is self funded with the possibility of funding from the Poliomyelitis Research Foundation
# Currently any licensing and usage of this software is governed under the regulations of the afore mentioned parties

#Author:	Keaghan Brown (3687524) - MSc Bioinformatics Candidate (3687524@myuwc.ac.za)
#Author:	Ruben Cloete (Supervisor) - Lecturer at South African National Bioinformatics Institute (ruben@sanbi.ac.za)

# This script completes the Second Objective for The Automated Computational Workflow to Prioritize Potential Resistance Variants Identified in HIV Integrase Subtype C and CRF02_AG
import os
import tkinter as tk 
from tkinter.filedialog import askopenfilename,  askdirectory
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Polypeptide import PPBuilder
import ntpath
import warnings
from Bio import BiopythonWarning
from get_raw_distances import *

class MutantAnalysis:
    def __init__(self):
        pymol.cmd.reinitialize()
        self.path = '/home/rotan/Desktop/Masters/mutant_pdbs'
        os.chdir(self.path)
        self.mutant_models = []
        for file in os.listdir(self.path):
            self.mutant_models.append(file)

    
    def residue_bonds(self):
        for model in self.mutant_models:
            parser = PDBParser(PERMISSIVE=1)
            structure_id = ntpath.basename(model).split('.')[0]
            filename = ntpath.basename(model)
            structure = parser.get_structure(structure_id, filename)
            ppb=PPBuilder()
            ppb.build_peptides(structure)
            for pp in ppb.build_peptides(structure):
                    sequence = pp.get_sequence()
            pymol.cmd.reinitialize()
            pymol.cmd.load(model)
            for i in cmd.get_chains(model):
                mut_resi_pos =  model.split(':')[1] 
                mut_resi_pos =  mut_resi_pos.split('.')[0] 
                mut_resi_pos =  mut_resi_pos[1:(len(mut_resi_pos)-1)]
                cmd. select('mutant_pos_'  + str(mut_resi_pos) ,  'resi ' +  str(mut_resi_pos))
                cmd.select('chain' + str(i),  'chain  ' + str(i))
                cmd.distance('h_bonds' + str(i),  "chain" + str(i), 'mutant_pos_' + str(mut_resi_pos), "4",  mode="2")
                D = get_raw_distances('h_bonds' + str(i))
                print(model)
                print(i)
                print("number of h-bonds:", len(D))
                

p = MutantAnalysis()
p.residue_bonds()
