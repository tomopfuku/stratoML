import tree_reader
import node
import sys

def get_tip_indices(tree):
    index = 0
    tipdic = {}
    for n in tree.iternodes():
        if n.istip:
            tipdic[n.label] = index
            index+=1
    return tipdic

def get_empty_bs(tree,typ="list"):
    if typ=="str":
        empty_bs = ""
        for n in tree.iternodes():
            if n.istip:
                empty_bs+="0"
    elif typ=="list":
        empty_bs = []
        for n in tree.iternodes():
            if n.istip:
                empty_bs.append("0")

    return empty_bs


def make_bipart(tree,node,tipdic):
    if node == tree:
        print("cannot make bipartition of the root")
        sys.exit()

    bs = get_empty_bs(tree,"list")
    for n in node.iternodes():
        if n.istip:
            nind = tipdic[n.label]
            bs[nind] = "1"
    return "".join(bs)


def decompose_tree(tree,tipdic):
    bp = []
    for n in tree.iternodes():
        if n == tree or len(n.children) == 0:
            continue

        bs = make_bipart(tree,n,tipdic)
        bp.append(bs)
    
    return frozenset(bp)


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("usage: "+sys.argv[0]+ " <tree 1 (newick)> <tree 2 (newick)>")
        sys.exit()
    
    nwk = open(sys.argv[1],"r").readline().strip()
    tree = tree_reader.read_tree_string(nwk)
    tipdic = get_tip_indices(tree)
    print(tree.get_newick_repr())
    print("TREE1:",tipdic)
    bp = decompose_tree(tree,tipdic)
    print("TREE1:",bp)

    nwk = open(sys.argv[2],"r").readline().strip()
    tree1 = tree_reader.read_tree_string(nwk)
    tipdic1 = get_tip_indices(tree)

    print("\n\nTREE2:",tipdic)
    bp1 = decompose_tree(tree1,tipdic1)

    print("TREE2:",bp1)

    z = bp.intersection(bp1)
    print("\noverlap:",z)

    y = bp.symmetric_difference(bp1)
    print("\nsymm. difference:",y)
    print("rf distance:",len(y))
    diff_res = abs(len(bp) - len(bp1))#/2.
    print("one tree is",diff_res,"biparts more resolved than the other")