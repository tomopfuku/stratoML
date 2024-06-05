import node,tree_reader,tree_utils
import stratlike
import sys


def reroot_tree(oldroot,newroot,bifroot = True):
    print(oldroot.get_newick_repr())
    if oldroot.istip == False: # need to adjust the root
        if len(oldroot.children) != 2:
            print("if root is hypothetical node, it needs to be bifurcating")
            print(newroot.get_newick_repr(),newroot.label())
            sys.exit()
        newpar = [ch for ch in oldroot.children if ch.subtree == newroot.subtree][0]
        newch = [ch for ch in oldroot.children if ch.subtree != newroot.subtree][0]
        newpar.parent = None
        newpar.children.append(newch)
        newch.parent = newpar
        oldroot.children = []
        print(newpar,newpar.label)
    else:
        newroot_node = node.Node()
    curnode = newroot.parent
    newroot.parent = oldroot
    oldroot.children.append(newroot)
    nextnode = curnode.parent
    curnode.parent = oldroot
    curnode.children.remove(newroot)
    oldroot.children.append(curnode)
    for _ in range(1000):
        lastnode = curnode
        curnode = nextnode
        nextnode = curnode.parent
        curnode.parent = lastnode
        lastnode.children.append(curnode)
        curnode.children.remove(lastnode)
        if nextnode == None:
            break

    print(oldroot.get_newick_repr())

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("usage:",sys.argv[0]," <tree> <stratigraphic ranges> <outgroup>")
        sys.exit()

    nwk = open(sys.argv[1],"r").readline().strip()
    tree = tree_reader.read_tree_string(nwk)
    tree_utils.map_strat_to_tree(tree,sys.argv[2])
    nodes = [n for n in tree.iternodes()]
    #print([n.label for n in nodes])
    #reroot_tree(tree,nodes[5])
    count = 0
    for n in tree.iternodes():
        n.index = count
        count+=1

    tree_utils.map_strat_to_tree(tree,sys.argv[2])
    print("Name,Code,Start,End,FAD,LAD,2.5HPD,97.5HPD,Parent")
    for n in tree.inorder():
        spls = [n.label,n.index,n.lower,n.upper,n.strat[0],n.strat[1],n.strat[0],n.strat[0]]
        if n != tree:
            spls.append(n.parent.index)
        else:
            spls.append("NA")
        spls = [str(i) for i in spls]
        print(",".join(spls))
