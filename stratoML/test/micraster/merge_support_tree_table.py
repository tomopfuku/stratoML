import sys

if len(sys.argv) != 3:
    print("usage: "+sys.argv[0]+" <support_table> <tree_table>")
    sys.exit()


suppfl = open(sys.argv[1],"r")
supp = {}
for line in suppfl:
    spls=line.strip().split()
    tax = spls[0]
    currsupp = spls[-1]
    if currsupp != "1.0":
        currsupp = str(round(float(currsupp),2))

    supp[tax] = currsupp

suppfl.close()

treefl = open(sys.argv[2],"r")
first = treefl.readline()
print(first.strip()+",support")
for line in treefl:
    spls=line.strip().split(",")
    tax = spls[0]
    currsupp = None
    try:
        currsupp = supp[tax]
    except:
        currsupp = ""
    spls.append(currsupp)
    print(",".join(spls))



