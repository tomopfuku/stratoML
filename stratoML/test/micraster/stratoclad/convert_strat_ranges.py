import sys


fl = open(sys.argv[1],"r")
h = fl.readline()
states = "0123456789ABCDEFGHIJKLMNOPQRSTUV"
statedic = {}
for i,j in enumerate(states):
    statedic[i]=j

ranges = {}
minl = 10000000
for line in fl:
    spls=line.strip().split()
    sp = spls[0]
    fa = int(spls[1])
    la = int(spls[2])
    if fa < minl:
        minl = fa

    ranges[sp] = (fa,la)    


for sp in ranges:
    newfa = ranges[sp][0] - minl
    newla = ranges[sp][1] - minl
    rangestr = "{"
    for i in range(newfa,newla+1):
        rangestr += statedic[i]
    rangestr += "}"
    print(sp,rangestr)
        

