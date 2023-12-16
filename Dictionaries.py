print ('imput password:')
pw = input('')
d = {}
print (pw)

for word in pw.split(' '):

    # Count the occurrences of the word in the original string
    count = pw.split().count(word)
    
    # Assign the count to the dictionary
    d[word] = count
    
for word, count in dict.items(d):
    print (word, count)
    
    
    
    
#     ## OR??? ##
    
# from collections import Counter
# for k,v in Counter(pw.split()).items(): print (k,v)