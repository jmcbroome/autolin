import bte
import sys
t = bte.MATree(sys.argv[1])
anns = t.dump_annotations()
for ann, node in anns.items():
    for s in t.get_leaves_ids(node):
        print(ann + "\t" + s)