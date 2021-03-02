# Generates random nucleotide sequence

import random

dna = ['A','C','G','T']
size = int(input(How many nucleotides would you like the random sequence to contain?\n"))
i = 0
seq = ""
while i < size:
    seq += random.choice(dna)
    i += 1
fileName = "randomDNAsequence" + str(size) + ".txt"
outputFile = open(fileName, "w")
outputFile.write("Randomly-generated sequence:\n")
outputFile.write(seq)
outputFile.close()
print ("Successfully saved generated sequence to \'" + fileName + "\'\n")

.............................................
OUTPUT
    
...........................................                 
