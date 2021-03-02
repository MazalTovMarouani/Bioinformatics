# Dictionary/Genome

#Problem 1
def readDataFromFile(fileName):# To get the file name. Return a dictionary. The keys are names, and the values are the associated read sequences.
    seq=open(fileName,'r')# Pointer to read the file 
    sequences=dict()
    for line in seq:# To read line by line 
        sequences[line[0:line.find(" ")]]=line[line.find(" ")+1:].replace('\n',"")#create an item in the dict--key is the name of the sequance,value is the sequence without the spaces
    return sequences # To return the dictionary

#Problem 2
def meanLength(fileName):# To get the file name returns the mean length of the reads 
    Sum=0.0
    i=0.0
    for key in seq:# To run on the dictionary that readDataFromFile returned-on each sequence
        Sum+=len(seq[key])# To sum the sequences lengths 
        i+=1# To sum the number of sequences
    return float(Sum/i)# To return the average

#Problem 3
def getOverlap(left, right):# To take 2 string arguments, left and right, each containing a read sequence. To return the overlapping sequence.
    overlap=""
    for i in range (0,len(left)):# To run on the left sequence
        if (left[i:]==right[:len(left)-i]):# To compare the largest possible overlap- get the smaller in each iteration 
            overlap=left[i:]# When we find the overlap
            break
    return overlap

#Problem 4
def getAllOverlaps(reads):# To compute overlaps between all pairs of reads in both left-right and right-left orientations.
    comp={}# Dictionary of dictionaries 
    for i in range(1,len(reads)+1):# To run on the number of sequences
        comp[str(i)]={}# To insert the name
        for j in range (1,len(reads)+1):# To run on the number of sequences , for the sequence that mentioned above,compute overlaps with all the others
            if (i!=j):
                comp[str(i)][str(j)]=len(getOverlap(reads[str(i)],reads[str(j)]))
    return comp

#Problem 5
def prettyPrint(overlaps):
    print "   ",
    for i in range(1, len(overlaps)+1):# To print the first raw of the sequences names
        print "% 3d" %i,
    print 
    for i in range(1,len(overlaps)+1): 
        print "% 3d" %i, # To print the relevant name 
        for j in range (1,len(overlaps)+1): 
            if (i==j):# Because there is no compare with the sequence itself 
                print "% 3s" % '-',
            if (i!=j):
                print "% 3d" %int (overlaps[str(i)][str(j)]),# To print the relevent overlap
        print

#Problem 6
def findFirstRead(overlaps):
    count=0
    for i in range(1,len(overlaps)+1):# To run on the keys in the total dictionary - the sequences names (left)
        for j in range (1,len(overlaps)+1):# To run on the keys in each inner dictionary- the sequences names(rigth)
            if (i!=j):
                if (overlaps[str(j)][str(i)]<2): # To check if the overlap is less than 2
                    count+=1
        if (count==len(overlaps)-1):# If all the overlaps of that sequences is less than 2
            return str(i)# To return the sequences name (left)
        else:
            count=0

#Problem 7
def findKeyForLargestValue(d):
    largest=0
    for key in d:# For each key in a given dicionary
        if (d[key]>largest):# To find the largest value
            largest= d[key]
            l=key
    return l,largest


def findOrder(name, overlaps):# To return a list with the order of the reads 
    nextName,largest=findKeyForLargestValue(overlaps[str(name)])# nextName=key of the lagest value. largest=largest value
    if (largest<2):# If we get the last read 
        return [name]
    else:
        return [name]+findOrder(str(nextName),overlaps)# To return recursively , and add to the list 

#Problem 8
def assembleGenome(readOrder, reads, overlaps):
    i=1# If it's the first read 
    j=0# To run on the reads number- allowing access to the previous read
    genome=""
    for name in readOrder:
        if (i==1):# If it's the first read 
            genome=reads[str(name)]
            i=0
        else:
            genome+=reads[name][overlaps[readOrder[j]][name]:]# To add to the genome the read without the overlap with the previous read 
            j+=1
    return genome
        

#main    
FileName=raw_input("Please enter a name of file (with .txt) ")# Gets the file name 
seq=readDataFromFile(FileName)# To get the file name. Return a dictionary. The keys are names, and the values are the associated read sequences.
print "The raw fragments:"
for key in seq:
    print key,":",seq[key]
Sum=meanLength(FileName)# To get the file name returns the mean length of the reads.
print "The mean length:"
print "%.2f"%Sum
comp=getAllOverlaps(seq)# Computes overlaps between all pairs of reads in both left-right and right-left orientations.
print "The overlaps matrix:"
prettyPrint(comp)# Prints the overlaps matrix 
s=findFirstRead(comp)# Returns the name of the first read
a=findOrder(str(s),comp)# Returns a list with the order of the reads
print "The order of the fragments:"
print a
genome=assembleGenome(a,seq,comp)# Components the genome
print "The final joined sequence:"
print genome


##OUTPUT
##Please enter a name of file (with .txt) genome assembly.txt
##The raw fragments:
##1 : GGCTCCCCACGGGGTACCCATAACTTGACAGTAGATCTCGTCCAGACCCCTAGC
##3 : GTCTTCAGTAGAAAATTGTTTTTTTCTTCCAAGAGGTCGGAGTCGTGAACACATCAGT
##2 : CTTTACCCGGAAGAGCGGGACGCTGCCCTGCGCGATTCCAGGCTCCCCACGGG
##5 : CGATTCCAGGCTCCCCACGGGGTACCCATAACTTGACAGTAGATCTC
##4 : TGCGAGGGAAGTGAAGTATTTGACCCTTTACCCGGAAGAGCG
##6 : TGACAGTAGATCTCGTCCAGACCCCTAGCTGGTACGTCTTCAGTAGAAAATTGTTTTTTTCTTCCAAGAGGTCGGAGT
##The mean length:
##55.33
##The overlaps matrix:
##      1   2   3   4   5   6
##  1   -   1   0   0   1  29
##  2  13   -   1   0  21   0
##  3   0   0   -   1   0   1
##  4   1  17   1   -   2   0
##  5  39   1   0   0   -  14
##  6   0   0  43   1   0   -
##The order of the fragments:
##['4', '2', '5', '1', '6', '3']
##The final joined sequence:
##TGCGAGGGAAGTGAAGTATTTGACCCTTTACCCGGAAGAGCGGGACGCTGCCCTGCGCGATTCCAGGCTCCCCACGGGGTACCCATAACTTGACAGTAGATCTCGTCCAGACCCCTAGCTGGTACGTCTTCAGTAGAAAATTGTTTTTTTCTTCCAAGAGGTCGGAGTCGTGAACACATCAGT
