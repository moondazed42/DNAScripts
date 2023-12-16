bpcomp = {'A':'T', 'T':'A', 'C':'G', 'G':'C', 'N':'N'}

while True:
    print("Enter Sequence:")
    call = input()
    rc = []

    for nuc in call:
        if nuc in bpcomp:
            rc.insert(0, bpcomp[nuc])
        else:
            print("Major Groove: Stop right there that's not a valid DNA sequence!")
            print ("Minor Groove: Major I found the culprit!--->",nuc)
            print ("Major Groove: Ah I see, well our friend here better fix that.")
            break
    else:
        # This block will execute if the for loop completes without encountering a 'break'
        print("Reverse Complement:", ''.join(rc))
        break  # Exit the while loop if the input was valid
