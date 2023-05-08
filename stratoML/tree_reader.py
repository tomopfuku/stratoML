import sys
import node

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
        newnode.lower = lower
        newnode.upper = upper
        newnode.length = lower - upper
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
    
    return tree

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("usage: "+ sys.argv[0]+ " <tree table file>")
        sys.exit()
    
    tree = read_tree_table(sys.argv[1])
    print(tree)
    print([n.label for n in tree.iternodes()])
    print(tree.get_newick_repr(True))