from Bio import Phylo
from io import StringIO
import matplotlib.pyplot as plt


tree1_1 = "(((1 :0.1, 2 :0.2) :0.12, 3 :0.3) :0.12, 4 :0.4);"
tree1_2 = "(((1 :0.1, 3 :0.2) :0.12, 2 :0.3) :0.12, 4 :0.4);"
tree1_3 = "(((1 :0.1, 2 :0.2) :0.12, 3 :0.3) :0.12, 4 :0.4);"

tree2 = "(((1 :0.1, 2 :0.2) :0.12, 3 :0.3) :0.12, 4 :0.4);"


tree1_1 = Phylo.read(StringIO(tree1_1), "newick")
tree1_2 = Phylo.read(StringIO(tree1_2), "newick")
tree1_3 = Phylo.read(StringIO(tree1_3), "newick")

tree2 = Phylo.read(StringIO(tree2), "newick")

# print(Phylo.draw_ascii(tree1_1))
# print(Phylo.draw_ascii(tree1_2))
# print(Phylo.draw_ascii(tree1_3))
# print(Phylo.draw_ascii(tree2))

# help(Phylo.draw)
Phylo.draw(tree1_1, axis=(False, ), savefig=('tree1_1.pdf', ))
Phylo.draw(tree1_2, axis=(False, ), savefig=('tree1_2.pdf', ))