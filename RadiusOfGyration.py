# Radius of Gyration

from Bio.PDB import *
import math
import glob
import warnings
warnings.filterwarnings("ignore")

HydrophilicResidues=["ASP","GLU","GLY", "HIS" ,"LYSV" ,"ASN" ,"GLN", "ARG", "SER", "THR", "TYR"]
HydrophobicResidues=["ALA","CYS","PHE","ILE","LEU","MET","PRO","VAL","TRP"]

class protein:
    ca=()
    re=()

    def __init__(self, FileName):
        parser = PDBParser()
        self.structure = parser.get_structure(FileName,FileName)

    def masscenter(self):
        sumX=0
        sumY=0
        sumZ=0
        for i in range(0,len(self.ca),3):
            sumX=sumX+self.ca[i]
            sumY = sumY + self.ca[i+1]
            sumZ = sumZ + self.ca[i+2]

        masscenter=[sumX/(len(self.ca)/3),sumY/(len(self.ca)/3),sumZ/(len(self.ca)/3)]
        return masscenter

    def Rg(self):
        MC=self.masscenter()
        #print("mass center:", MC)
        d=0
        for i in range(0,len(self.ca),3):
            d=d+(((self.ca[i]-MC[0])**2)+((self.ca[i+1]-MC[1])**2)+((self.ca[i+2]-MC[2])**2))**0.5
        d=d/len(self.re)
        #print("Rg",d,len(myProtein.ca)),
        return d

    def Rg_phobic(self):
        MC = self.masscenter()
        d = 0
        count2=0
        for i in range(0, len(self.re)):
            if self.re[i] in HydrophobicResidues:
                count2 = count2 + 1
                d = d + (((self.ca[i*3] - MC[0]) ** 2) + ((self.ca[i*3 + 1] - MC[1]) ** 2) + ((self.ca[i*3 + 2] - MC[2]) ** 2))**0.5
        d = d / count2
        #print("p", d)
        return d
 def Rg_philic(self):
        MC = self.masscenter()
        d = 0
        count = 0
        count2 = 0
        for i in range(0, len(self.re)):
            if self.re[i] in HydrophilicResidues:
                count2 = count2 + 1
                d = d + (((self.ca[i*3] - MC[0]) ** 2) + ((self.ca[i*3 + 1] - MC[1]) ** 2) + ((self.ca[i*3 + 2] - MC[2]) ** 2))**0.5
        d = d / count2
        #print("h", d)
        return d



for pdbfile in glob.glob('*.pdb'):
    myProtein=protein(pdbfile)
   #print("the name of file is:",str(pdbfile))
    for model in myProtein.structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    if atom.name=="CA":
                        x,y,z=atom.get_coord()
                        myProtein.ca=myProtein.ca+(x,y,z)
                        name=residue.get_resname()
                        myProtein.re=myProtein.re+(name,)
    print("Rg: %8.2f  Rg phobic: %8.2f  Rg philic:  %8.2f  length: %8.2f"% (myProtein.Rg(), myProtein.Rg_phobic(), myProtein.Rg_philic(),len(myProtein.re)))
    #myProtein.Rg()
    #myProtein.Rg_phobic()
    #myProtein.Rg_philic()
#print (myProtein.ca)
#print (len(myProtein.ca))
#print (len(myProtein.ca))
#print (myProtein.re)
#print (len(myProtein.re))
