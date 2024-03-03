import sys,re


fl = open(sys.argv[1],"r")

names = {}
reading_names = False
tree_block = False
for line in fl:
    if line.strip() == "":
        continue
    if line.strip()=="BEGIN TREES;":
        tree_block = True
        continue
    if line.strip()=="TRANSLATE":
        reading_names = True
        continue
    if reading_names and line.strip() != ";":
        spls = line.strip().split()
        names[spls[0]] = spls[1].replace(",","")
    elif reading_names and line.strip() == ";":
        reading_names = False
        continue
    if tree_block == True and reading_names == False and line.strip()!=";":
        nwk = line.strip().split()[-1]
        for code in reversed(range(len(names))):
            ccode = str(code+1)
            nwk = re.sub(ccode,names[ccode],nwk)
        print(nwk)
    elif tree_block == True and reading_names == False and line.strip()==";":
        break

