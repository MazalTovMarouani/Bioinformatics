# ORF

def ORF (dna,pixORF): #get dna string, and a number that complete to the length of the dna , return ORFs
    flag=0# updated when we found ORF
    i=0
    x=0
    for i in range(0,len(dna),3):#find the first "ATG"
        wordstart= dna[i:i+3]#read every time 3 letters 
        if (wordstart=="ATG"):# if word="ATG"
            break
    while ( i<len(dna)-299):# until we find "ATG" and the length will be less then 300 that then it dosen't fit to the demands
        x=i # x= to the location of "ATG"
        for x in range( i,len(dna),3):#find the first stop codon after the "ATG"
            word = dna[x:x+3]#read every time 3 letters 
            if ((word=="TAG")or(word=="TGA")or (word=="TAA")):#if word= any stop codon 
                if((x+3)-i>300):# if the length between the start and the end that found is more then 300
                   print "start:",i+pixORF,"end:",x+pixORF+2,"length:", (x+3)-i,dna[i:i+6],dna[x-3:x+3]
                   flag=1
                   print dna[i:x+3]# print the ORF
                   ProteinTranslation (dna[i:x+3])#send the ORF to tranlate to protein
                   break # exsit from the stop codon searching 
                else:
                    break# exsit from the stop codon searching 
        i=x # updated the start location to start from the end codon that found 
        for n in range(i,len(dna) ,3):# searching for new start codon 
            wordstart= dna[n:n+3]#read every time 3 nucleotides
            if (wordstart=="ATG"):
                break
        i=n # updated the start location to the new location of "ATG"
    if (flag==0):# if we did'nt find any ORF
        print "ERROR, ORF,doesn't exist"

def ProteinTranslation (ORF):#translate a given ORF to preducted protein
    # dictionary that containing  the codes for every amino acid 
    code = {"TTT":"F", "TTC":"F", "TTA":"L", "TTG":"L","TCT":"S", "TCC":"s", "TCA":"S", "TCG":"S","TAT":"Y", "TAC":"Y", "TAA":"STOP", "TAG":"STOP","TGT":"C", "TGC":"C", "TGA":"STOP", "TGG":"W",
    "CTT":"L", "CTC":"L", "CTA":"L", "CTG":"L", "CCT":"P", "CCC":"P", "CCA":"P", "CCG":"P","CAT":"H", "CAC":"H", "CAA":"Q", "CAG":"Q","CGT":"R", "CGC":"R", "CGA":"R", "CGG":"R",
    "ATT":"I", "ATC":"I", "ATA":"I", "ATG":"M","ACT":"T", "ACC":"T", "ACA":"T", "ACG":"T" ,"AAT":"N", "AAC":"N", "AAA":"K", "AAG":"K","AGT":"S", "AGC":"S", "AGA":"R", "AGG":"R",
    "GTT":"V", "GTC":"V", "GTA":"V", "GTG":"V","GCT":"A", "GCC":"A", "GCA":"A", "GCG":"A","GAT":"D", "GAC":"D", "GAA":"E", "GAG":"E","GGT":"G", "GGC":"G", "GGA":"G", "GGG":"G"} 
    protein=""
    for i in range(0,len(ORF),3):#for every 3 nucleotides, add to 'protein' the suitable amino acid
        word=ORF[i:i+3]#read every time 3 nucleotides
        protein=protein+code[word]# add to 'protein' the suitable amino acid
    print "Predicted protein:", protein

def Reverse (dna):# get a string and return the inverse complementary string 
    compliteDna=""
    ReverseDna=""
    x=0
    code={"T":"A","A":"T","G":"C","C":"G"} # dictionary that containing the complementary nucleotide for all the nucleotides
    for i in dna:# for each nucleotide in 'dna' 
        compliteDna= compliteDna+code[i]# add the complementary nucleotide 
    while (x<len( compliteDna)): # for each nucleotide in  the complementary sequence 
        ReverseDna=ReverseDna+compliteDna[len(compliteDna)-1-x]# write form the end to the start 
        x=x+1
    return ReverseDna

      

FileName=raw_input("please enter name of file (with .txt) ")#get the file name
DNA=file(FileName,"r").read()# reads the file to DNA
DNA=DNA.replace(DNA[DNA.find(">"):DNA.find("\n")],"")#remove the first line 
DNA=DNA.replace("\n","")#remove all the \n 
DNA=DNA.upper()#convert to capital letters
# print which ORF is that and send to the func suitability
print "ORF 1+ : "
ORF(DNA,1)
print "ORF 2+ :"
ORF(DNA[1:],2)
print "ORF 3+ :"
ORF(DNA[2:],3)
ReverseDNA=Reverse(DNA)# return the complementary nucleotide
print "ORF 1- :"
ORF(ReverseDNA,1)
print "ORF 2- :"
ORF(ReverseDNA[1:],2)
print "ORF 3- :"
ORF(ReverseDNA[2:],3)
