# GC content

import random
def bases_freq(dna):#calculates each nucleotid frequency and return as dictionary
    lenr=0
    nuc=dict()
    nuc['A']=0
    nuc['T']=0
    nuc['C']=0
    nuc['G']=0
    nuc['GC']=0
    seq=dna.read(100)# read only 100 chars each time
    while(seq):
        lenr+=len(seq)-seq.count("\n")# calculate the length
        nuc['A']+=(seq.count('A'))
        nuc['T']+=(seq.count('T'))
        nuc['C']+=(seq.count('C'))
        nuc['G']+=(seq.count('G'))
        nuc['GC']+=(seq.count('C')+(seq.count('G')))
        seq=dna.read(100)
        seq.upper() 
    nuc['A']= nuc['A']/float (lenr) #divide with the length 
    nuc['T']=nuc['T']/float (lenr)
    nuc['C']=nuc['C']/float (lenr)
    nuc['G']=nuc['G']/float (lenr)
    nuc['GC']=nuc['GC']/float (lenr)
    dna.seek(0,0)#returns to the file's start
    return nuc

def slidingwindowplot(dna,window_length):#the func get file and window length,returns two lists
    #one with the window's value, the second with the gc content
    x=[]
    y=[]
    i=0
    seq=dna.read(window_length)#read only the window length each time 
    while (seq):
       y.append((seq.count('C')+seq.count('G'))/float(window_length))#add to the list the frequency of gc 
       x.append (i+window_length)
       i=i+window_length # promote i 
       seq=dna.read(window_length)
    dna.seek(0,0)#returns to the file's start
    return x,y
def representation(word, dna_seq_file):#calculates whether a particular string is more or less than we would expect to see by chance
    gc=0
    length=0
    seq=dna.read(10000)
    while (seq):
        gc+=seq.count(word)
        length+=len(seq)-seq.count("\n")-seq.count("N")
        seq=dna.read(10000)
    dna_seq_file.seek(0,0)
    x=bases_freq(dna_seq_file)
    freq=float(gc)/length
    a=x['G']*x['C']
    return freq/a
    
rand=open("random.txt","w")
for i in range(20000):
    rand.write(random.choice("ATGC"))
FileName=raw_input("please insert the file's name(with .txt) ")
DNA= open( FileName,"r")
NumNuc=bases_freq(DNA)
print FileName
print " a   t     c     g    gc"
print"%.1f"%(NumNuc['A']*100),"%.1f"%(NumNuc['T']*100),"%.1f"%(NumNuc['C']*100),"%.1f"%(NumNuc['G']*100),"%.1f"%(NumNuc['GC']*100)
Tmp=open("tmp.txt","w")# open new file 
for line in DNA:# copy only the DNA Sequence without opening line and without spaces
    x=line.find(">")
    if (x==-1):
        line.replace("N","")
        Tmp.write("".join(line.split("\n")))
    
Tmp.close()
dna=open("tmp.txt","r")
x,y=slidingwindowplot(dna,input("please enter a window_length "))
ex=open("GCexel.txt","w")
ex.write("x: ")
for i in x:# write x's values to the file 
    ex.write(str(i))
    ex.write(" ")
ex.write("\ny: ")
for i in y:# write y's values to the file 
    ex.write(str(i))
    ex.write(" ")
Tmp1=open("tmp.txt","r")
GCcontent=representation(raw_input("plese enter a 'word'  in capital letters "), Tmp1) 
print GCcontent   
DNA.close()
