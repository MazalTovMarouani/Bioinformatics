#Mazale-Tov MAROUANI ID 323315721

#function to edit the text
def editText (seq):
    seq = seq.strip(seq.split("\n")[0])    #delete the begin of the fasta formt
    seq = seq.replace("\n","")              #delete the whitwspaces
    return seq

#function to build the matrix
def matrix (seq1, seq2):
    a = []                        #Defining the matrix
    for i in range(len(seq2)+2):        #build the rows
        a.append([])                    #build the columns
        for j in range(len(seq1)+2):
            if ((i == 0) and (j == 0)):
                a[i].append(" ")      #fill the upper left corner 
            elif (i == 0):           #fill the first row
                if (j == 1):
                    a[i].append("0")      
                else:
                    a[i].append(seq1[j-2])
            elif (j == 0):           #fill the first column
                if (i == 1):
                    a[i].append("0")
                else:
                    a[i].append(seq2[i-2])
            else:                     #fill the body of the matrix
                a[i].append(" ")
    return a

#function to fill the alignment and pointers matrix
def fillMatrix (aliMat, poiMat, match, mis, gap):
    for i in range(1,len(seq2)+2):     #fill the first column
        aliMat[i][1] = gap*(i-1)      #fill the alignment matrix
        poiMat[i][1]  = 'V'           #fill the pointers matrix
    for j in range(1,len(seq1)+2):      #fill the first row
        aliMat[1][j] = gap*(j-1)
        poiMat[1][j]  = 'H'
    poiMat[1][1]  = '0'               #fill the upper left corner 
    for i in range(2,len(seq2)+2):
        for j in range(2,len(seq1)+2):
            if (aliMat[i][0] == aliMat[0][j]):       #if the letters are the same
                d = aliMat[i-1][j-1] + match       #the score in case the current cell is computed from the upper left cell
            else:                                 #if the letters are not the same
                 d = aliMat[i-1][j-1] + mis
            v = aliMat[i-1][j] + gap            #the score in case the current cell is computed from the upper cell
            h = aliMat[i][j-1] + gap          #the score in case the current cell is computed from the left cell
            cell = max(d, v, h)              #check what is the max score
            aliMat[i][j] = cell
            if (cell == d):     #fill the pointers matrix
                poiMat[i][j] = 'D'
            elif (cell == h):
                poiMat[i][j] = 'H'
            elif (cell == v):
                poiMat[i][j] = 'V'
    mat = [aliMat,poiMat]    #to return the two matrix
    return mat


#function to print a matrix
def printMat(matrix, rows, columns):
    for i in range(rows):
        print
        for j in range(columns):
            print "% 3s" % matrix[i][j],
    print "\n"


#recursive function that returns the path 
def path (matrix, strPath, i, j):
    if (matrix[i][j] == "0"):       #Stop conditions
        return strPath
    else:
        if (matrix[i][j] == 'D'):   #in case the current cell is computed from the upper left cell
            strPath = matrix[i][j] + path (matrix, strPath, i-1, j-1)  #call the recursive function, and send the upper left cell
        elif (matrix[i][j] == 'H'):     #in case the current cell is computed from the upper cell
            strPath = matrix[i][j] + path (matrix, strPath, i, j-1)
        else:                         #in case the current cell is computed from the left cell
            strPath = matrix[i][j] + path (matrix, strPath, i-1, j)
    return strPath

#function to find and print the actual pairwise alignment
def alignment (seq1, seq2, strPath):
    for i in range(len(strPath)):      
        if (strPath[i] == "H"):   #in case there is a gap at seq2
            seq2 = seq2[:i] + "-" + seq2[i:]  #enter "-" to seq2 at the place of the gap
        elif (strPath[i] == "V"):    #in case there is a gap at seq1
            seq1 = seq1[:i] + "-" + seq1[i:]  #enter "-" to seq1 at the place of the gap
    print seq1
    print seq2
    return
        
#main program
nameFile = raw_input("please enter file name: ")
file = open(nameFile, 'r')
dna = file.read()      #read the dna sequences to a string variable
file.close()
dna = dna.split(">")     #Divide the string into two sequences
seq1 = dna[1]
seq2 = dna[2]
seq1 = editText(seq1)   #edit the sequnce
seq2 = editText(seq2)   #edit the sequnce
match = int(raw_input("please enter match score: "))   #the score of match
mis = int(raw_input("please enter mismatch score: "))   #the score of mismatch
gap = int(raw_input("please enter gap score: "))     #the score of gap

aliMat = matrix(seq1, seq2)     #build the alignment matrix  
poiMat = matrix(seq1, seq2)     #build the pointers matrix
mat = fillMatrix (aliMat, poiMat, match, mis, gap)   #fill the alignment and the pointers matrixes
aliMat = mat[0]
poiMat = mat[1]
print "Alignment Matrix:"
printMat(aliMat, len(seq2)+2, len(seq1)+2)    #print the alignment matrix
print "Pointers Matrix:"
printMat(poiMat, len(seq2)+2, len(seq1)+2)   #print the pointers matrix
strPath = ""
strPath = path (poiMat, strPath, (len(seq2)+1), (len(seq1)+1)) #use function to find the path
print "Path:\n"
print strPath, "\n"   #print the path
print "Pairwise Alignment:\n"
alignment (seq1, seq2, strPath)   #find and print the actual pairwise alignment



