import sys
import numpy as np
import smaps

def recode_poly_traits(traits,ss):
    retraits = {}
    for n in traits:
        curtr = traits[n]
        new = []
        for i,tr in enumerate(curtr):
            if tr[0] == -9:
                new.append(-9)
                continue
            curss = ss[i]
            if curss < 2:
                new.append(1)
                continue
            smap = smaps.get_smap(curss)
            retr = -9
            for ii,row in enumerate(smap):
                #print(list(row))
                allpres = True
                for j in tr:
                    if row[j] == 0:
                        allpres = False
                if allpres == True and sum(row) == len(tr):
                    retr = ii
            new.append(retr)
            #print(tr,retr)
        retraits[n] = new
        #print(n,new)
    return retraits

def read_fasta(flnm):
    fl = open(flnm,"r", encoding="utf-8")
    traits = {}
    allstates = {} 
    for line in fl:
        if len(line.strip()) == 0:
            continue
        if line.strip()[0] == ">":
            curlab = line.strip().replace(">","")
        else:
            linels = line.strip().split()
            seq = [[0]]
            for i,trait in enumerate(linels):
                tt = trait.strip().split("|")
                states = []
                for ttt in tt:
                    if ttt == "?" or ttt == "-":
                        states.append(-9)
                    else:
                        try:
                            allstates[i].append(int(ttt))
                        except:
                            allstates[i] = []
                            allstates[i].append(int(ttt))

                        states.append(int(ttt))
                seq.append(states)
            traits[curlab] = seq
            #print(curlab,seq)
    state_spaces = [2]
    for i in range(len(allstates)):
        curst = set(allstates[i])
        state_spaces.append(max(curst)+1)
    return traits, np.array(state_spaces,dtype=int)


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("usage: "+ sys.argv[0]+ " <fasta file>")
        sys.exit()
    
    seq = read_fasta(sys.argv[1])
    print(seq)
