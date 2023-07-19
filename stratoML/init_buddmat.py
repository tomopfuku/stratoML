import smaps
import numpy as np
#import split_mat
import sys

def init_buddmat():
    all_comp = []
       
    """
    for nstate in range(2,4):
        cur = smaps.get_smap(nstate)
        for i0 in range(len(cur)):
            par = cur[i0]
            print(i0,list(par))
    """
                   
    print("import numpy as np\n\n")
    decl = "cdef long[:,:,:] "
    start = 2
    stop = 5
    for i in range(start,stop):
        decl += f"buddmat{i}"
        if i < stop-1:
            decl+=","

    #print(decl+"\n")

    for nstate in range(start,stop):
        statecomp = []
        cur = smaps.get_smap(nstate)
        for i0 in range(len(cur)):
            if i0 == 0:
                statecomp.append([])
                continue
            par = cur[i0]
            chrg = []
            if sum(par) == 1:
                chrg.append(i0)
                statecomp.append(chrg)
                continue
            for i1 in range(len(cur)):
                if i1 == 0:
                    continue
                ch1 = cur[i1]
                
                if sum(ch1) <= sum(par):
                    nogood=False
                    for j in range(len(par)):
                        if par[j] == 0 and ch1[j] == 1:
                            nogood = True
                            break
                    if nogood==False:
                        chrg.append(i1)
            statecomp.append(chrg)
        all_comp.append(statecomp)

    
    allt = []
    for i in range(len(all_comp)):
        nstate = i + 2
        curbudd = all_comp[i]
        nmax = len(curbudd[-1])
        t = np.zeros((len(curbudd),nmax),dtype=int)
        for j in range(len(curbudd)):
            curst = curbudd[j]
            for k in range(len(curst)):
                t[j][k] = curbudd[j][k]
        print("buddmat"+str(nstate)+" = np.array(\n"+str([list(m) for m in list(t)])+",dtype=int)\n")
        allt.append(t)
    #print("np.array(\n"+str(allt))

    print("def get_buddmat(nstates):")
    #print("def get_spltmat(int nstates):")
    for i in range(start,stop):
        print(f"    if nstates == {i}:\n        return buddmat{i}")



init_buddmat()
