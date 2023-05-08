import sys
import numpy as np

def read_fasta(flnm):
    fl = open(flnm,"r")
    traits = {}
    allstates = {} 
    for line in fl:
        if len(line.strip()) == 0:
            continue
        if line[0] == ">":
            curlab = line.strip().replace(">","")
        else:
            linels = line.strip().split()
            seq = []
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
    state_spaces = []
    for i in range(len(allstates)):
        curst = set(allstates[i])
        state_spaces.append(len(curst))
    return traits, np.array(state_spaces,dtype=int)


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("usage: "+ sys.argv[0]+ " <fasta file>")
        sys.exit()
    
    seq = read_fasta(sys.argv[1])
    print(seq)