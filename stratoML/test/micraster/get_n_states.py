import sys

fl = open(sys.argv[1],"r")
tax = ""
tax_nstate = {}
for line in fl:
    if line[0]==">":
        tax = line.strip().replace(">","")
    else:
        chars = line.strip().split()
        ntrait = len(chars)
        nstate = 0
        for i in chars:
            spls = i.strip().split("|")
            nstate += len(spls)
        tax_nstate[tax] = str(nstate-ntrait)
        #print(tax,str(nstate-ntrait))

fl = open(sys.argv[2],"r")
h = fl.readline()
pars = {}
codes = {}
for line in fl:
    spls=line.strip().split(",")
    parent = spls[-1]
    code = spls[1]
    tax = spls[0]
    pars[tax] = parent
    codes[code] = tax

print("taxon nstates par_nstates")
for tax in tax_nstate:
    parent = codes[pars[tax]]
    if parent == "":
        par_nstates = "NA"
    else:
        par_nstates = tax_nstate[parent]
    print(tax,tax_nstate[tax],par_nstates)