import sys
import node
import numpy as np

def read_tree_table(flnm):
    fl = open(flnm,"r")
    h = fl.readline()
    codes = {}
    tree = None
    for line in fl:
        newnode = node.Node()
        spls = line.strip().split(",")
        label = spls[0]
        code  = spls[1]
        lower = float(spls[2])
        upper = float(spls[3])
        parent = spls[4]
        newnode.label = label
        if label != "":
            newnode.istip = True
        newnode.lower = lower + 0.01
        if upper > 0.01:
            newnode.upper = upper - 0.01
        else:
            newnode.upper = upper 
        #newnode.length = newnode.lower - newnode.upper
        newnode.strat = np.array([lower,upper],dtype = np.double)
        if parent == "NA":
            tree = newnode

        codes[code] = [parent,newnode]

    for ncode in codes:
        par = codes[ncode][0]
        if par == "NA":
            continue
        curnode = codes[ncode][1]
        parnode = codes[par][1]
        parnode.add_child(curnode)
        
    for n in tree.iternodes(0):
        if n == tree:
            continue

        if n.istip == False:
            oldest_ch = max([c.lower for c in n.children])
            n.upper = oldest_ch
        else:
            if len(n.children) > 0:
                youngest_ch = min([c.lower for c in n.children])
                if n.upper > youngest_ch:
                    n.upper = youngest_ch

        if n.lower < n.parent.upper:
            n.lower = n.parent.upper

        n.length = n.lower - n.upper
        #print(n.label,n.lower,n.upper,n.strat[0],n.strat[1])
    return tree

def init_heights_strat(tree,fixed_root=False):
    for i in tree.iternodes(order=1):
        if i == tree and fixed_root == True:
            continue
        elif i.istip and i.num_occurrences != 0: #Set heights for tips
            if i.lower != 0.0:
                i.height = i.strat[1] - 0.01
            else:
                i.height = i.strat[1]
        elif i.istip == False:# and i.label == "":
            o = []
            start = False
            for j in i.children:
                if j.num_occurrences != 0:
                    if start == False:
                        o = [j.upper,j.lower]
                        #print o
                        start = True
                    else:
                        o = o+[j.upper,j.lower]
            if len(o) > 0:
                i.height = max(o+[j.height for j in i.children]) + 0.01
            else:
                i.height = max([j.height for j in i.children])+0.01
    for i in tree.iternodes():
        if i == tree:
            continue
        i.length = i.parent.height-i.height



if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("usage: "+ sys.argv[0]+ " <tree table file>")
        sys.exit()
    
    tree = read_tree_table(sys.argv[1])
    print(tree)
    print([n.label for n in tree.iternodes()])
    print(tree.get_newick_repr(True))
