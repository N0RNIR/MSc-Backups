import os;
def foldx_analysis(self, pdb_path):
    self.pdb_path = pdb_path
    os.chdir(self.pdb_path )
    for file in os.listdir(self.pdb_path):
            if file.endswith('.pdb'):
    foldx = '/home/rotan/Desktop/FoldX/foldx --command=Stability  --pdb=s_test.pdb'
    os.system(foldx)