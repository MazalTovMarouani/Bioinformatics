#To find the most frequent k-mers in a string.

f=open("dna.txt")
kmers=dict() 
k=(input("please enter a value: k="))
r=f.read()
length=((len(r))-k+1)
for i in range (length): 
    if(r[0:k] in kmers):#if this word already exists
        kmers[r[0:k]]=kmers[r[0:k]]+1#increments the value of this key 
    else:
        kmers[r[0:k]]=1#new word:value=1
    r=(r[1:len(r)])#prepare the next check
frequent=max(kmers.values())
for key in kmers: 
    if(kmers[key]==frequent):
       print (key)
       print(kmers[key])#to print the words of kmer with the highest frequency
f.close()
