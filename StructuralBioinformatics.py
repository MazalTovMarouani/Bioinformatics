#Structural bionformatics

from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio import SeqIO
import numpy as np

class seq(ProteinAnalysis):
    charges = {'D':-1, 'E':-1, 'H':1, 'C':-1, 'Y':-1, 'K':1, 'R':1, 'N-ter':1, 'C-ter':-1}
    pKa={'D':4.05, 'E':4.45, 'H':5.98, 'C':9.0, 'Y':10.0, 'K':10.0, 'R':12.0, 'N-ter':7.5, 'C-ter':3.55}

    def proteinCharge(self,pH):
        '''to calculate the charge of an entire protein at a specific pH'''
        charge=0
        pka=0
        #charges of the edges:
        charge=charge+(1/(1+10**(pH-self.pKa.get('N-ter'))))#N terminal
        charge=charge+(-1)/(1+10**(self.pKa.get('C-ter')-pH))#C terminal
        
        for res in self.sequence: 
            c=self.charges.get(res)
            if (c!=None):
                if (c==1):#basic aa
                    pka=self.pKa.get(res)
                    charge=charge+(1)/(1+10**(pH-pka))
                if (c==-1):#acid aa
                    pka=self.pKa.get(res)
                    charge=charge+(-1)/(1+10**(pka-pH))
            else:
                charge+=0 
        return charge
    
    def pI(self):
         '''to calculate the pI of a given protein'''
         pI=0
         charge=500
         pH_list=np.arange(0,14,0.1) #list which contain pHs  from 0 to 14 in intervals of 0.1
         for i in pH_list: #for each pH in the list
             c=self.proteinCharge(i)
             if(abs(c)<abs(charge)):  
                 charge=c
                 pI=i
         return pI
        

     
    def Molecular_Weight(self):
        '''return the molecular mass of a given protein'''
        return self.molecular_weight()

    def ResidueNumber(self):
        count=0
        res=self.count_amino_acids()#res is dictionnary wich gives us the number of each aa that there is in the protein seq
        for i in res.values(): #"for" that pass on each value in the dic res
            count=count+i
        return count
    

# main:

File = open("protSeq.txt", "r")
f=SeqIO.parse(File, "fasta")
Seq = SeqIO.to_dict(f)
File.close()

for v in Seq.values():
    print(v.id)
    print(v.seq)
    p=""
    p=str(v.seq)
    print("pI: "+str(seq(p).pI()))
    print("Number of Amino acids: "+str(seq(p).ResidueNumber()))
    print("Protein molecular mass: "+str(seq(p).Molecular_Weight())+" Da")+("\n")

'''output:
Cow_ubi
MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG
pI: 6.8
Number of Amino acids: 76
Protein molecular mass: 8564.7357 Da

Human_ubi
MQIFVKTRKGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG
pI: 9.1
Number of Amino acids: 76
Protein molecular mass: 8634.8322 Da

Zebra-fish_ubi
MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKDEEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG
pI: 5.3
Number of Amino acids: 78
Protein molecular mass: 8808.9371 Da

Chimpanzee_ubi
MQIFVKTLETGKTITLEVEPSDTIENVKAKIQDKEGIPPEEDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG
pI: 5.1
Number of Amino acids: 79
Protein molecular mass: 8952.0777 Da
'''
