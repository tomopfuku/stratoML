import numpy as np


def all_bs(nstates):
    ls = []
    for i in range(1<<nstates):
        bs = bin(i)[2:]
        bs = '0'*(nstates-len(bs))+bs
        lbs = [int(j) for j in list(bs)]
        lbs.reverse()
        ls.append(lbs)
    return ls


if __name__ == "__main__":
    print("import numpy as np\n\n")
    decl = "cdef long[:,:] "
    start = 2
    stop = 8
    for i in range(start,stop):
        decl += f"smap{i}"
        if i < stop-1:
            decl+=","

    print(decl+"\n")

    for i in range(start,stop):
        bs = all_bs(i)
        print(f"smap{i} = np.array([")
        for row in bs:
            print(str(row)+",")
        print("],dtype=int)\n")

    print("def get_smap(int nstates):")
    for i in range(start,stop):
        print(f"    if nstates == {i}:\n        return smap{i}")

