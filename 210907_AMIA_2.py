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
import itertools
from tabulate import tabulate


class MutantAnalysis:
    def __init__(self):
        pymol.cmd.reinitialize()
       
    
    def residue_bonds(self, mutant, model):
        self.mutant_path = mutant
        self.model_path = model
        os.chdir(self.mutant_path)
        for mutant in os.listdir(self.mutant_path):
            if mutant .endswith('.pdb'):
                table_dict = {"Model Contacts" : [],  "Mutant Contacts" : [],  "Contacts Gained/ Lost":  []}
                mutant_contacts = {}
                model_contacts = {}
                print(mutant.split('.')[0].replace(':', '_'))
                parser = PDBParser(PERMISSIVE=1)
                structure_id = ntpath.basename(mutant).split('.')[0]
                filename = ntpath.basename(mutant)
                structure = parser.get_structure(structure_id, filename)
                ppb=PPBuilder()
                ppb.build_peptides(structure)
                cmd.reinitialize()
                cmd.load(mutant)
                mut_resi_pos =  mutant.split(':')[1] 
                mut_resi_pos =  mut_resi_pos.split('.')[0] 
                mut_resi_pos =  mut_resi_pos[1:(len(mut_resi_pos)-1)]
                cmd. select('mutant_pos_'  + str(mut_resi_pos) ,  'model ' + str(mutant.split('.')[0].replace(':', '_') )+ '& resi ' +  str(mut_resi_pos))
                chains = cmd.get_chains(mutant)
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
                        cmd.distance('mutant_h-bond' +  str(i)+str(j),  'mutant_pos_'  + str(mut_resi_pos), 'residue' + str(i)+str(j), "3.5", mode="2")
                        D = get_raw_distances('mutant_h-bond' +  str(i)+str(j))
                        mutant_contacts[str(i)+ ','+str(j)] = (len(D))
                cmd.reinitialize()
                cmd.load(self.model_path)
                cmd. select('model_pos_'  + str(mut_resi_pos) ,  'model ' + str(ntpath.basename(self.model_path).split('.')[0])+ '& resi ' +  str(mut_resi_pos))
                for k in range(len(largest_seq)):
                    for l in chains:
                        cmd.select('model_residue' + str(k)+str(l), 'resi ' +  str(k) + '& chain ' + str(l))
                        cmd.distance('model_h-bond' +  str(k)+str(l),  'model_pos_'  + str(mut_resi_pos), 'model_residue' + str(k)+str(l), "3.5", mode="2")
                        D2 = get_raw_distances('model_h-bond' +  str(k)+str(l))
                        model_contacts[str(k)+ ','+ str(l)] = (len(D2))
                for (key, key2)  in  zip(mutant_contacts, model_contacts):
                    contact_diff = mutant_contacts[key] - model_contacts[key2]
                    if contact_diff < 0 or contact_diff > 0:
                        table_dict.setdefault('Contacts Gained/ Lost', []).append(key + ': ' + str(contact_diff))
                for  contact in model_contacts:
                    if model_contacts[contact] > 0:
                        table_dict.setdefault('Model Contacts', []).append(contact + ': ' + str(model_contacts[contact] ))
                for  contact2 in mutant_contacts:
                    if mutant_contacts[contact2] > 0:
                        table_dict.setdefault('Mutant Contacts', []).append(contact2 + ': ' + str(mutant_contacts[contact2] ))
                data_file = open('/home/rotan/Desktop/Objective2/' + str(mutant.split('.')[0].replace(':', '_')) + '_contacts.txt', 'w')
                data_file.write(tabulate(table_dict, headers="keys"))
                data_file.close()
             
    def foldx_analysis(self, pdb_path):
        self.pdb_path = pdb_path
        os.chdir(self.pdb_path)
        for file in os.listdir(self.pdb_path):
                if file.endswith('.pdb'):
                    foldx_repair = '/home/rotan/Desktop/FoldX/foldx --command=RepairPDB --pdb=' + file
                    foldx_stability = '/home/rotan/Desktop/FoldX/foldx --command=Stability  --pdb=' + file  + ' --output-dir=/home/rotan/Desktop/Objective2/Obj2.P2_ouput/' 
                    os.system(foldx_repair)
                    os.system(foldx_stability)

p = MutantAnalysis()
p.foldx_analysis('/home/rotan/Desktop/Masters/mutant_pdbs')

