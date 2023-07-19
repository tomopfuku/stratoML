import smaps
import numpy as np
#import split_mat

def init_splitmat():
    all_comp = []
    """
    for nstate in range(2,4):
        cur = smaps.get_smap(nstate)
        for i0 in range(len(cur)):
            par = cur[i0]
            print(i0,list(par))
    """
    print("from numpy import array\n\n")
    decl = "cdef long[:,:,:] "
    start = 2
    stop = 5
    for i in range(start,stop):
        decl += f"spltmat{i}"
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
                chrg.append([i0,i0])
                statecomp.append(chrg)
                continue
            for i1 in range(len(cur)):
                if i1 == 0:
                    continue
                ch1 = cur[i1]
                for i2 in range(len(cur)):
                    if i2 == 0:
                        continue
                    ch2 = cur[i2]
                    if sum(ch1)+sum(ch2) == sum(par) or sum(ch1)+sum(ch2) < sum(par)*2 or sum(ch1)+sum(ch2) == sum(par)*2:
                        numover = 0
                        nogood=False
                        for j in range(len(par)):
                            if par[j] == 1 and ch1[j] == 0 and ch2[j] == 0:
                                nogood = True
                            elif par[j] == 0 and (ch1[j] == 1 or ch2[j] == 1):
                                nogood = True
                                break
                        if nogood==False:
                            chrg.append([i1,i2])
            statecomp.append(chrg)
        all_comp.append(statecomp)

    allt = []
    for i in range(len(all_comp)):
        nstate = i + 2
        curbif = all_comp[i]
        nmax = len(curbif[-1])
        t = np.zeros((len(curbif),nmax,2),dtype=int)
        for j in range(len(curbif)):
            curst = curbif[j]
            for k in range(len(curst)):
                t[j][k] = curbif[j][k]
        print("spltmat"+str(nstate)+" = array(\n"+str(list(t))+")\n")
        allt.append(t)
    #print("np.array(\n"+str(allt))

    print("def get_spltmat(nstates):")
    #print("def get_spltmat(int nstates):")
    for i in range(start,stop):
        print(f"    if nstates == {i}:\n        return spltmat{i}")



init_splitmat()
