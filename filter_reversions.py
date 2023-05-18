import sys
sys.path.append("~/bin:")
import bte
treename = sys.argv[1]
t = bte.MATree(treename)
print("Tree loaded.")
if sys.argv[4] != None:
    thold = int(sys.argv[4])
else:
    print("Using default threshold (2 reversions on a branch)")
    thold = 2
to_prune = set()
with open(sys.argv[2]) as inf:
    first = True
    for entry in inf:
        if first:
            first = False
            continue
        node, lc, mc, md, reversions = entry.strip().split()
        if int(reversions) >= int():
            if lc == '1':
                to_prune.add(node)
            else:
                descendents = t.get_leaves_ids(node)
                for d in descendents:
                    to_prune.add(node)
print(f"Identified {len(to_prune)} samples to remove.")
to_keep = [s for s in t.get_leaves_ids() if s not in to_prune]
sub = t.subtree(to_keep)
sub.save_pb(sys.argv[3])
print("Filtering complete.")