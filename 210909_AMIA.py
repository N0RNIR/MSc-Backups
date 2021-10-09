# This script is developed for the fufuillment for Masters at the South African National Bioinformatics Institute at the University of the Western Cape
# Project Title: "Automated computational workflow to prioritize potential resistance variants identified in HIV Integrase Subtype C and CRF02_AG" 
# Currently the project is self funded with the possibility of funding from the Poliomyelitis Research Foundation
# Currently any licensing and usage of this software is governed under the regulations of the afore mentioned parties

#Author:	Keaghan Brown (3687524) - MSc Bioinformatics Candidate (3687524@myuwc.ac.za)
#Author:	Ruben Cloete (Supervisor) - Lecturer at South African National Bioinformatics Institute (ruben@sanbi.ac.za)

# This script completes the first objective for The Automated Computational Workflow to Prioritize Potential Resistance Variants Identified in HIV Integrase Subtype C and CRF02_AG

import pymol
from Bio.PDB.PDBParser import PDBParser 
from Bio.PDB.Polypeptide import PPBuilder
import ntpath
import warnings
from Bio import BiopythonWarning
warnings.simplefilter('ignore', BiopythonWarning)

class Mutations:
    def __init__(self):
        pymol.cmd.reinitialize()
        self.three_letter ={'V': 'VAL', 'I': 'ILE', 'L': 'LEU', 'E': 'GLU', 'Q': 'GLN', 'D': 'ASP', 'N': 'ASN', 'H': 'HIS', 'W': 'TRP', 'F': 'PHE', 'Y': 'TYR', 'R': 'ARG', 'K': 'LYS', 'S': 'SER', 'T': 'THR', 'M': 'MET', 'A': 'ALA', 'G': 'GLY', 'P': 'PRO', 'C': 'CYS'}
        self.drug_mutation_dictionary = {} 
    
    def process_mutant_list(self,  mutant_list):
        self.mutant_list = mutant_list
        my_file = open(self.mutant_list)
        file_data = my_file.read()
        mutation_data = file_data.split('\n')                                                                                                                                      
        while '' in mutation_data:                                                                                                                                                       
            mutation_data.remove('')
        mutation_list_counter = 0                                                                                                                                                      
        for j in mutation_data:                                                                                                                                                          
            mutation_list_counter += 1                                                                                                                                               
            if j == '###':                                                                                                                                                                          
                drug_name = mutation_list_counter                                                                                                                          
            elif j == "####":                                                                                                                                                                   
                first_mutation = mutation_list_counter                                                                                                                       
            elif j == "##":                                                                                                                                                                       
                last_mutation = mutation_list_counter - 1                                                                                                                  
                self.drug_mutation_dictionary[mutation_data[drug_name]] = (str(first_mutation) + ','  + str(last_mutation))       
        for key in self.drug_mutation_dictionary:                                                                                                                                  
            start_stop_pos = self.drug_mutation_dictionary[key].split(',')
            self.drug_mutation_dictionary[key] = mutation_data[int(start_stop_pos [0]) : int(start_stop_pos [1])] 
    
    def single_mutations(self, model):
        self.model=model
        structure_id = self.model
        file_name = self.model
        parser = PDBParser(PERMISSIVE = 1)
        structure = parser.get_structure(structure_id, file_name)
        ppbuilder = PPBuilder()
        start_pos = ppbuilder.build_peptides(structure)[0][0].get_id()[1]
        for key in self.drug_mutation_dictionary:
            for mutation in  self.drug_mutation_dictionary[key]:
                pymol.cmd.reinitialize()
                initial_residue = mutation[0]
                mutated_residue = mutation[len(mutation)-1]
                residue_pos = mutation[1:len(mutation)-1] 
                residue_pos = start_pos + int(residue_pos)-1
                print(mutation)
                cmd.load(self.model)
                for i in cmd.get_chains(ntpath.basename(self.model).split('.')[0]):
                    cmd.wizard('mutagenesis')
                    cmd.refresh_wizard()
                    cmd.get_wizard().do_select(i + '/' + str(residue_pos) + '/')
                    cmd.get_wizard().set_mode(self.three_letter[mutated_residue])
                    cmd.get_wizard().apply()
                    cmd.set_wizard()
                cmd.save('/home/rotan/Desktop/Masters/mutant_pdbs/'  + str(key) + mutation + '.pdb') 
                
    def multi_mutations(self, model):
        self.model=model
        structure_id = self.model
        file_name = self.model
        parser = PDBParser(PERMISSIVE = 1)
        structure = parser.get_structure(structure_id, file_name)
        ppbuilder = PPBuilder()
        start_pos = ppbuilder.build_peptides(structure)[0][0].get_id()[1]
        cmd.load(self.model)
        for key in self.drug_mutation_dictionary:
            for mutation in  self.drug_mutation_dictionary[key]:
                initial_residue = mutation[0]
                mutated_residue = mutation[len(mutation)-1]
                residue_pos = mutation[1:len(mutation)-1] 
                residue_pos = start_pos + int(residue_pos)-1
                print(mutation)
                for i in cmd.get_chains(ntpath.basename(self.model).split('.')[0]):
                    cmd.wizard('mutagenesis')
                    cmd.refresh_wizard()
                    cmd.get_wizard().do_select(i + '/' + str(residue_pos) + '/')
                    cmd.get_wizard().set_mode(self.three_letter[mutated_residue])
                    cmd.get_wizard().apply()
                    cmd.select('mutation', 'resi ' + str(residue_pos) + ' around 3')
                    cmd.clean('mutation')
                    cmd.set_wizard()
        cmd.save('/home/rotan/Desktop/Masters/mutant_pdbs/Multiple_Mutants.pdb')

p = Mutations()
p.process_mutant_list("/home/rotan/Desktop/mut_lst.txt")
p.multi_mutations('/home/rotan/Desktop/210814_S_4_mer_Model.pdb')