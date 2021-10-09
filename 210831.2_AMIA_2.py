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
    
    def residue_bonds(self, path):
        self.path = path
        os.chdir(self.path)
        for model in os.listdir(self.path):
            print(model)
            parser = PDBParser(PERMISSIVE=1)
            structure_id = ntpath.basename(model).split('.')[0]
            filename = ntpath.basename(model)
            structure = parser.get_structure(structure_id, filename)
            ppb=PPBuilder()
            ppb.build_peptides(structure)
            pymol.cmd.reinitialize()
            pymol.cmd.load(model)
            mut_resi_pos =  model.split(':')[1] 
            mut_resi_pos =  mut_resi_pos.split('.')[0] 
            mut_resi_pos =  mut_resi_pos[1:(len(mut_resi_pos)-1)]
            cmd. select('mutant_pos_'  + str(mut_resi_pos) ,  'resi ' +  str(mut_resi_pos))
            chains = cmd.get_chains(model)
            sequences = []
            for pp in ppb.build_peptides(structure):
                sequence= pp.get_sequence()
                sequences.append(sequence)
            largest_seq =  ''
            for i in  sequences:
                if len(i) > len(largest_seq):
                    largest_seq = i
            for i in range(len(largest_seq)):
                for j in chains:
                    cmd.select('residue' + str(i)+str(j), 'resi ' +  str(i) + '& chain ' + str(j))
                    cmd.distance('h-bond' +  str(i)+str(j),  'mutant_pos_'  + str(mut_resi_pos), 'residue' + str(i)+str(j), "3.5", mode="2")
                    D = get_raw_distances('h-bond' +  str(i)+str(j))
                    print('h-bond' +  str(i)+str(j))
                    print("number of h-bonds:", len(D))


p = MutantAnalysis()
p.residue_bonds('/home/rotan/Desktop/Masters/mutant_pdbs')
