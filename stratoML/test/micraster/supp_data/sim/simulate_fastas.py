import sys
sys.path.append("../../")
import tree_reader,tree_utils
import node
import qmat
import buddmat as bm
import stratlike
import numpy as np
from simulate_mfc import *
#from random import choice
from random import choices
from random import random
import smaps


if len(sys.argv) != 6:
    print("usage:" + sys.argv[0] + " <tree> <stratigraphic data> <num_traits> <max_num_states> <out_DIR>")
    sys.exit()

nwk = open(sys.argv[1],"r").readline()
tree = tree_reader.read_tree_string(nwk)
tree_utils.map_strat_to_tree(tree,sys.argv[2]) 
num_traits = int(sys.argv[3])
ss = int(sys.argv[4])
for rep in range(100):
    rand_m = 0.02 + ( random() * (0.2 - 0.01) )
    #rand_m = 0.05
    #rand_f = 0.01 + ( random() * (0.1 - 0.01) )
    rand_f = rand_m
    qmats = qmat.Qmat(rand_m,rand_f)
    stratlike.calibrate_brlens_strat(tree,0.3)
  
    tree_utils.sort_children_by_age(tree) 
    tree_utils.init_budd_marginals(tree,num_traits,ss)

    sim_traits_across_tree(tree,qmats,num_traits,ss)
    trait_dict = get_trait_dict(tree,ss)
    outfl = open(sys.argv[5]+"/"+str(rep)+"_"+str(round(rand_m,3))+"_"+str(round(rand_f,3))+".fa","w")

    for i in trait_dict:
        outfl.write(">"+i+"\n")
        outfl.write(" ".join(trait_dict[i])+"\n")
    outfl.close()
