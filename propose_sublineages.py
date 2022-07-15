import bte
import sys
import argparse
def simple_node_distance(t, nid, pnid):
    td = 0
    for anc in t.rsearch(nid,True):
        if anc.id == pnid:
            return td
        td += len(anc.mutations)
    return td

def evaluate_candidate(t, a, nid, pgp_d, ignore = set()):
    """Evaluate a candidate branch as a putative sublineage.

    Args:
        t (MATree): The tree.
        a (str): The parent lineage annotation node.
        nid (str): The node id of the candidate branch.
    """
    leaves = [l for l in t.get_leaves(nid) if l.id not in ignore]
    if len(leaves) == 0:
        return 0
    candidate_to_parent = simple_node_distance(t, nid, a)
    if candidate_to_parent == 0:
        return 0
    total_distances = 0
    for l in leaves:
        dist = simple_node_distance(t, l.id, nid)
        total_distances += dist
    if total_distances == 0:
        return 0
    mean_distances = total_distances/len(leaves)
    candidate_value = len(leaves) * (candidate_to_parent / (mean_distances + candidate_to_parent))
    mean_distances_parent = (candidate_to_parent*len(leaves) + total_distances)/len(leaves)
    parent_value = len(leaves) * (pgp_d / (mean_distances_parent + pgp_d))
    return candidate_value - parent_value

def get_plin_distance(t,nid):
    td = 0
    for n in t.rsearch(nid,True):
        td += len(n.mutations)
        try:
            if len(n.annotations) > 0 or n.is_root():
                return td
        except:
            continue
    return td

def evaluate_lineage(t, anid, ignore = set(), floor = 0, maxpath = 100):
    """Evaluate every descendent branch of lineage a to propose new sublineages.

    Args:
        t (MATree): The tree.
        a (str): The lineage annotation node to check.
    """
    parent_to_grandparent = min(get_plin_distance(t,anid), maxpath)
    candidates = t.depth_first_expansion(anid)
    good_candidates = []
    for c in candidates:
        if not c.is_leaf():
            cscore = evaluate_candidate(t, anid, c.id, parent_to_grandparent, ignore) - floor
            if cscore > 0:
                good_candidates.append((cscore,c))
    if len(good_candidates) == 0:
        return (0,None)
    return max(good_candidates, key=lambda x: x[0])

def argparser():
    parser = argparse.ArgumentParser(description="Propose sublineages for existing lineages based on relative representation concept.")
    parser.add_argument("-i", "--input", required=True, help='Path to protobuf to annotate.')
    parser.add_argument("-c", "--clear", action='store_true', help='Clear all current annotations and apply a level of serial annotations to start with.')
    parser.add_argument("-r", "--recursive", action='store_true', help='Recursively add additional sublineages to proposed lineages.')
    parser.add_argument("-o", "--output", help='Path to output protobuf, if desired.',default=None)
    parser.add_argument("-d", "--dump", help="Print proposed sublineages to a table.",default=None)
    parser.add_argument("-l", "--labels", help="Print samples and their associated lowest level lineages to a table.",default=None)
    parser.add_argument("-f", "--floor", help="Gain of a proposed and current lineage label must be more than this much. Default 0",type=float,default=0)
    parser.add_argument("-m", "--maxpath", help="Set a maximum path length value when computing sublineage viability. Reduce to allow clades descended from long branches to be further subdivided. Default 100",type=int,default=100)
    args = parser.parse_args()
    return args

def main():
    args = argparser()
    t = bte.MATree(args.input)
    if args.dump != None:
        dumpf = open(args.dump,'w+')
    if args.clear:
        t.apply_annotations({node.id:[] for node in t.depth_first_expansion()})
    annotes = t.dump_annotations()
    if len(annotes) == 0:
        print("No lineages found in tree; starting from root.")
        annotes = {'L':t.root.id}
    # print("Tree contains {} annotated lineages initially.".format(len(annotes)),file=sys.stderr)
    #keep going until the length of the annotation dictionary doesn't change.
    if args.dump != None:
        print("parent\tparent_nid\tproposed_sublineage\tproposed_sublineage_nid\tproposed_sublineage_score",file=dumpf)
    outer_annotes = annotes
    while True:
        new_annotes = {}
        for ann,nid in outer_annotes.items():
            all_leaves = t.get_leaves_ids(nid)
            print("Considering descendents of node {} with annotation {}, with {} total leaves.".format(nid,ann,len(all_leaves)),file=sys.stderr)
            serial = 0
            labeled = set()
            while True:
                best_score, best_node = evaluate_lineage(t, nid, ignore = labeled, floor = args.floor, maxpath = args.maxpath)
                if best_score <= 0:
                    break
                new_annotes[ann + "." + str(serial)] = best_node.id
                if args.dump != None:
                    print("{}\t{}\t{}\t{}\t{}".format(ann,nid,ann + "." + str(serial),best_node.id,str(best_score+args.floor)),file=dumpf)
                for l in t.get_leaves_ids(best_node.id):
                    labeled.add(l)
                if len(labeled) >= len(all_leaves):
                    break
                serial += 1
                print(serial)
        if not args.recursive:
            break
        elif len(new_annotes) == 0:
            break
        else:
            annotes.update(new_annotes)
            outer_annotes = new_annotes
    print("After sublineage annotation, tree contains {} annotated lineages.".format(len(annotes)),file=sys.stderr)
    if args.output != None:
        t.apply_annotations({v:[k] for k,v in annotes.items()})
        t.save_pb(args.output)
    if args.dump != None:
        dumpf.close()
    if args.labels != None:
        labels = {}
        for lid in t.get_leaves_ids():
            for n in t.rsearch(lid,True):
                try:
                    if len(n.annotations) > 0:
                        labels[lid] = n.annotations[0]
                        break
                except IndexError:
                    continue
        with open(args.labels,'w+') as f:
            print("sample\tlineage",file=f)
            for k,v in labels.items():
                print("{}\t{}".format(k,v),file=f)

if __name__ == "__main__":
    main()