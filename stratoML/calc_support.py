import sys, math
import node,tree_reader
import biparts as bp

def check_all_bps(tree_bp,cur_bp):
    if len(cur_bp) != len(tree_bp):
        return False
    elif len(cur_bp.symmetric_difference(tree_bp)) > 0:
        return False
    elif len(cur_bp.symmetric_difference(tree_bp)) == 0:
        return True
    

if len(sys.argv) != 3:
    print("usage: "+ sys.argv[0] + " <tree to map to> <list of trees with AIC scores>")
    sys.exit()

maptree = tree_reader.read_tree_string(open(sys.argv[1],"r").readline())

fl = open(sys.argv[2],"r")
trees = []
for line in fl:
    spls = line.strip().split()
    aic = float(spls[0])
    nwk = spls[1]
    tree = tree_reader.read_tree_string(nwk)
    trees.append((aic,tree))

best = 100000000000.
besttree = None
for i in trees:
    if i[0] < best:
        best = i[0]
        besttree = i[1]


tipdic = bp.get_tip_indices(maptree)
curbp = bp.decompose_tree(maptree,tipdic)
origbp = bp.decompose_tree(maptree,tipdic)
bpsupp = {}
for i in curbp:
    bpsupp[i]=0

bpall = []
sumrellike = 0.
origAIC = None
for i in trees:
    aic = i[0]
    rellike = math.exp(0.5*(best-aic))
    sumrellike+=rellike
    #print(rellike)
    curtree = i[1]
    curbp = bp.decompose_tree(curtree,tipdic)
    same_tree = check_all_bps(origbp,curbp)
    if same_tree:
        origAIC = rellike
    bpall.append((rellike,curbp))

for best_bp in bpsupp:
    cumaic = 0.
    for i in bpall:
        rellike = i[0]
        decomp_tree = i[1]
        for curbp in decomp_tree:
            #print(best_bp,curbp)
            if best_bp == curbp:
                weight = rellike / sumrellike
                bpsupp[best_bp] += weight

#print(bpsupp)

tree_supp = origAIC / sumrellike
print("TREE_SUPPORT:",tree_supp)

for n in maptree.iternodes():
    if n == tree or len(n.children) == 0:
        continue

    bs = bp.make_bipart(besttree,n,tipdic)
    for i in bpsupp:
        if i == bs:
            print(n.label,bs,bpsupp[i])

