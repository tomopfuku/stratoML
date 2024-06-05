import tree_reader
import sys

fl = sys.argv[1]
tree = tree_reader.read_tree_table(fl)
print(tree.get_newick_repr()+";")
