import bte
import sys
import argparse

def process_mstr(mstr):
    """Read a mutation string and return the chromosome, location, reference, and alternate alleles.
    """
    if ":" in mstr:
        chro = mstr.split(":")[0]
        data = mstr.split(":")[1]
    else:
        chro = None
        data = mstr
    if data[0].isdigit():
        loc = int(data[:-1])
        ref = None
        alt = data[-1]
    else:
        loc = int(data[1:-1])
        ref = data[0]
        alt = data[-1]
    return chro, loc, ref, alt

def dists_to_root(tree, node, mutweights = {}):
    #nodes must be a dict that gets updated on each recursion
    #gives back a dict with all nodes and their respective dist from root
    #initalize this with our starting node at 0, its our "root" whether its the actual tree root or not
    nodes = {node.id:0}
    def recursive_dists_to_roots(node):
        if node.id == tree.root.id:
            nodes[node.id] = 0   #becomes 0 because it is the root
        for child in node.children:
            if (node.id == tree.root.id):
                dist = sum([mutweights.get((int(m[1:-1]),m[-1]),1) for m in child.mutations])
            else:
                dist = nodes[node.id] + len(child.mutations)    
            nodes[child.id] = dist
            recursive_dists_to_roots(child)
    recursive_dists_to_roots(node)
    return nodes

def get_node_length(node, mutweights = {}):
    tlen = 0
    for m in node.mutations:
        _, loc, _, alt = process_mstr(m)
        tlen += mutweights.get((loc,alt),1)
    return tlen

def get_sum_and_count(rbfs, ignore = set(), mutweights = {}):
    # node sum stored in first index and node count stored in second index of each dict entry
    sum_and_count_dict = {}
    leaf_count = 0
    for node in rbfs:
        if node.is_leaf():
            leaf_count += 1
            if node.id not in ignore:
                sum_and_count_dict[node.id] = (get_node_length(node,mutweights), 1)
        else:
            total_count = 0
            total_sum = 0
            for child in node.children:
                sumtc = sum_and_count_dict.get(child.id, None)
                if sumtc == None:
                    continue
                total_count += sumtc[1]
                total_sum += sumtc[0]
            if total_count > 0:
                #total path length is computed as the total path lengths to each child plus the length of the current node TIMES the number of samples.
                #this is because total path length is not the same as tree parsimony- some mutations are part of many sample paths
                #for a given sample to its parent, the total path length is just the number of mutations (as computed above)
                #but for an internal node with two leaf children's path length with respect to its parent, 
                #its equal to the sum of the two child's path lengths plus 2 times its mutations, since those mutations are shared among 2 samples
                #this logic applies as we move further up the tree.
                sum_and_count_dict[node.id] = (total_sum + get_node_length(node,mutweights) * total_count, total_count)
    return sum_and_count_dict, leaf_count #, leaves

def evaluate_candidate(a, nid, sum_and_counts, dist_to_root, minimum_size=0,minimum_distinction=0):
    """Evaluate a candidate branch as a putative sublineage.

    Args:
        t (MATree): The tree.   
        a (str): The parent lineage annotation node.
        nid (str): The node id of the candidate branch.
    """
    node_sum, node_count = sum_and_counts.get(nid,[0,0])
    if node_sum == 0 or node_count == 0:
        return 0
    candidate_to_parent = dist_to_root[nid] - dist_to_root[a]
    mean_distances = node_sum/node_count
    if (mean_distances + candidate_to_parent) == 0:   #avoid divide by 0
        candidate_value = 0
    else:
        candidate_value = max([(node_count-minimum_size),0]) * (max([candidate_to_parent-minimum_distinction,0])) / (mean_distances + candidate_to_parent)
    return candidate_value

def get_plin_distance(t,nid,mutweights = {}):
    td = 0
    for n in t.rsearch(nid,True):
        td += get_node_length(n, mutweights)
        if any([ann != "" for ann in n.annotations]):
            return td
    return td

def evaluate_lineage(t, dist_to_root, anid, candidates, sum_and_count, minimum_size = 0, minimum_distinction = 0, banned = set()):
    """Evaluate every descendent branch of lineage a to propose new sublineages.

    Args:
        t (MATree): The tree.
        a (str): The lineage annotation node to check.
    """
    good_candidates = []
    for c in candidates:
        if not c.is_leaf() and c.id not in banned:
            cscore = evaluate_candidate(anid, c.id, sum_and_count, dist_to_root,minimum_size,minimum_distinction)
            if cscore > 0:
                good_candidates.append((cscore,c))
    if len(good_candidates) == 0:
        return (0,None)
    return max(good_candidates, key=lambda x: x[0])

def get_outer_annotes(t, annotes):
    """Get all outer annotations in a tree.

    Args:
        t (MATree): The tree.
        annotes (dict): The annotation dictionary.
    """
    #find the outermost annotation nodes by looking first at all annotations, then excluding ones that are 
    #ancestors to some other one. There's also no point in checking ancestor lineages for any lineage that is not itself outermost
    #as they are not outermost by definition, so use a skip list.
    skip = set()
    outer_annotes = {k:v for k,v in annotes.items()}
    for lin, nid in annotes.items():
        if lin in skip:
            continue
        ancestors = t.rsearch(nid, False)
        for anc in ancestors:
            for ann in anc.annotations:
                if ann != "":
                    outer_annotes.pop(ann, None)
                    skip.add(ann)
    return outer_annotes

def parse_mutweights(mutweights_file):
    """Parse a mutation weight file.
    """
    mutweights = {}
    with open(mutweights_file) as f:
        for line in f:
            line = line.strip()
            if line == "":
                continue
            if line[0] == "#":
                #ignore comment lines
                continue
            parts = line.split()
            _, loc, _, alt = process_mstr(parts[0])
            mutweights[(loc, alt)] = float(parts[1])
    return mutweights

def argparser():
    parser = argparse.ArgumentParser(description="Propose sublineages for existing lineages based on relative representation concept.")
    parser.add_argument("-i", "--input", required=True, help='Path to protobuf to annotate.')
    parser.add_argument("-c", "--clear", action='store_true', help='Clear all current annotations and apply a level of serial annotations to start with.')
    parser.add_argument("-r", "--recursive", action='store_true', help='Recursively add additional sublineages to proposed lineages.')
    parser.add_argument("-o", "--output", help='Path to output protobuf, if desired.',default=None)
    parser.add_argument("-d", "--dump", help="Print proposed sublineages to a table.",default=None)
    parser.add_argument("-l", "--labels", help="Print samples and their associated lowest level lineages to a table.",default=None)
    parser.add_argument("-t", "--distinction", help="Require that lineage proposes have at least i mutations distinguishing them from the parent lineage or root.",type=int,default=1)
    parser.add_argument("-m", "--minsamples", help="Require that each lineage proposal describe at least f samples.", type=int, default=10)
    parser.add_argument("-w", "--mutweights", help="Path to an optional two column space-delimited containing mutations and weights to assign them.",default=None)
    parser.add_argument("-g", "--gene", help='Consider only mutations in the indicated gene. Requires that --gtf and --reference be set.', default=None)
    parser.add_argument("-s", "--missense", action='store_true', help="Consider only missense mutations. Requires that --gtf and --reference be set.")
    parser.add_argument("--gtf", help="Path to a gtf file to apply translation. Use with --reference.")
    parser.add_argument("--reference", help='Path to a reference fasta file to apply translation. Use with --gtf.')
    args = parser.parse_args()
    return args

def main():
    args = argparser()
    t = bte.MATree(args.input)
    if args.gtf != None and args.reference != None:
        print("Performing tree translation and removing mutations not included in selection.")
        print("Initial tree parsimony:", t.get_parsimony_score())
        translation = t.translate(fasta_file=args.reference,gtf_file=args.gtf)
        to_use = {}
        for nid, aav in translation.items():
            for_nid = []
            for aa in aav:
                if args.missense and aa.is_synonymous():
                    continue
                if args.gene != None:
                    if aa.gene != args.gene:
                        continue
                for_nid.append(aa.nuc)
            to_use[nid] = for_nid
        t.apply_mutations(to_use)
        print("Mutations filtered; new tree parsimony:", t.get_parsimony_score())
    mutweights = {}
    if args.mutweights != None:
        mutweights = parse_mutweights(args.mutweights)
        
    if args.dump != None:
        dumpf = open(args.dump,'w+')
    if args.clear:
        t.apply_annotations({node.id:[] for node in t.depth_first_expansion()})
    annotes = t.dump_annotations()
    original_annotations = set(annotes.keys())
    if len(annotes) == 0:
        print("No lineages found in tree; starting from root.")
        annotes = {'L':t.root.id}
    else:
        print("{} annotations found in the tree; identifying candidates for subdivision.".format(len(annotes)))
        annotes = get_outer_annotes(t, annotes)
        print("{} outer annotations found in the tree; identifying sublineages.".format(len(annotes)))
    print("Tree contains {} annotated lineages initially.".format(len(annotes)),file=sys.stderr)
    #keep going until the length of the annotation dictionary doesn't change.
    if args.dump != None:
        print("parent\tparent_nid\tproposed_sublineage\tproposed_sublineage_nid\tproposed_sublineage_score",file=dumpf)
    outer_annotes = annotes
    level = 1
    while True:
        print("Level: ",level)
        new_annotes = {}
        used_nodes = set()
        for ann,nid in outer_annotes.items():
            serial = 0
            labeled = set()
            rbfs = t.breadth_first_expansion(nid, True) #takes the name
            dist_root = dists_to_root(t, t.get_node(nid), mutweights) #needs the node object, not just the name
            while True:
                scdict, leaf_count = get_sum_and_count(rbfs, ignore = labeled, mutweights = mutweights)
                best_score, best_node = evaluate_lineage(t, dist_root, nid, rbfs, scdict, args.minsamples, args.distinction, used_nodes)
                if best_score <= 0:
                    break
                newname = ann + "." + str(serial)
                for anc in t.rsearch(best_node.id,True):
                    used_nodes.add(anc.id)
                new_annotes[newname] = best_node.id
                if args.dump != None:
                    print("ADDING LINEAGES FROM LEVEL {}".format(level))
                    print("{}\t{}\t{}\t{}\t{}".format(ann,nid,newname,best_node.id,str(best_score)),file=dumpf)
                for l in t.get_leaves_ids(best_node.id):
                    labeled.add(l)
                
                if len(labeled) >= leaf_count:
                    break
                serial += 1
                print("Annotated lineage", serial)
        if not args.recursive:
            annotes.update(new_annotes)
            break
        elif len(new_annotes) == 0:
            break
        else:
            annotes.update(new_annotes)
            outer_annotes = new_annotes
            level += 1
    print("After sublineage annotation, tree contains {} annotated lineages.".format(len(annotes)),file=sys.stderr)
    if args.output != None:
        annd = {}
        for k,v in annotes.items():
            if v not in annd:
                annd[v] = []
            if len(annd[v]) == 2:
                annd[v][1] = k
            else:
                annd[v].append(k)
        #reload the original tree to restore all edits.
        t = bte.MATree(args.input)
        t.apply_annotations(annd)
        t.save_pb(args.output)
    if args.dump != None:
        dumpf.close()
    if args.labels != None:
        if args.output == None:
            annd = {}
            for k,v in annotes.items():
                if v not in annd:
                    annd[v] = []
            if len(annd[v]) == 2:
                annd[v][1] = k
            else:
                annd[v].append(k)
            t.apply_annotations(annd)
        labels = {}
        for ann, nid in annotes.items():
            if ann not in original_annotations:
                for l in t.get_leaves_ids(nid):
                    labels[l] = ann
        with open(args.labels,'w+') as f:
            print("strain\tlineage",file=f)
            for k,v in labels.items():
                print("{}\t{}".format(k,v),file=f)

if __name__ == "__main__":
    main()