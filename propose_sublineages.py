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

def compute_mutation_weight(node,mutweights):
    if len(mutweights) == 0:
        return len(node.mutations)
    dist = 0
    for m in node.mutations:
        _, loc, _, alt = process_mstr(m)
        mweight = max([mutweights.get((loc,alt,None),0),mutweights.get((loc,alt,node.id),0)])
        dist += mweight
    return dist

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
                dist = compute_mutation_weight(child,mutweights)
            else:
                dist = nodes[node.id] + compute_mutation_weight(child,mutweights)  
            nodes[child.id] = dist
            recursive_dists_to_roots(child)
    recursive_dists_to_roots(node)
    return nodes

def get_sum_and_count(rbfs, ignore = set(), mutweights = {}):
    # node sum stored in first index and node count stored in second index of each dict entry
    sum_and_count_dict = {}
    leaf_count = 0
    for node in rbfs:
        if node.is_leaf():
            leaf_count += 1
            if node.id not in ignore:
                sum_and_count_dict[node.id] = (compute_mutation_weight(node,mutweights), 1)
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
                sum_and_count_dict[node.id] = (total_sum + compute_mutation_weight(node,mutweights) * total_count, total_count)
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
        # print("DEBUG: {} {} {} {}".format(node_count, max([(node_count-minimum_size+1),0]), candidate_to_parent, max([candidate_to_parent-minimum_distinction+1,0])))
        candidate_value = max([(node_count-minimum_size+1),0]) * (max([candidate_to_parent-minimum_distinction+1,0])) / (mean_distances + candidate_to_parent)
    return candidate_value

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

def get_skipset(t, annotes):
    """
    Return the set of nodes which are, or are ancestral to, existing lineages on the tree. 
    Used to build a banned node list to prevent retroactive definition of lineage parents. 
    """
    skip = set()
    for lin, nid in annotes.items():
        ancestors = t.rsearch(nid, True)
        for anc in ancestors:
            skip.add(anc.id)
    return skip

def get_outer_annotes(t, annotes):
    """Get all outer annotations (annotations which are terminal for at least one sample) in a tree.

    Args:
        t (MATree): The tree.
        annotes (dict): The annotation dictionary.
    """
    outer_annotes = {}
    for l in t.get_leaves():
        mann = l.most_recent_annotation()
        for a in mann:
            if a != None and a not in outer_annotes:
                outer_annotes[a] = annotes[a]
        if len(outer_annotes) == len(annotes):
            #all of them will be checked. No need to continue.
            break
    return outer_annotes

def parse_mutweights(mutweights_file):
    """Parse a mutation weight file. First column is the mutation, second column is the weight. Optionally set a third column to be the occurrence nodes at which the weight is used.
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
            if len(parts) == 2:
                mutweights[(loc, alt, parts[2])] = float(parts[1])
            else:
                mutweights[(loc, alt, None)] = float(parts[1])
    if len(mutweights) == 0:
        print("ERROR: Mutation weight file indicated found empty!")
        exit(1)
    return mutweights

def argparser():
    parser = argparse.ArgumentParser(description="Propose sublineages for existing lineages based on relative representation concept.")
    parser.add_argument("-i", "--input", required=True, help='Path to protobuf to annotate.')
    parser.add_argument("-c", "--clear", action='store_true', help='Clear all current annotations and apply a level of serial annotations to start with.')
    parser.add_argument("-r", "--recursive", action='store_true', help='Recursively add additional sublineages to proposed lineages.')
    parser.add_argument("-o", "--output", help='Path to output protobuf, if desired.',default=None)
    parser.add_argument("-d", "--dump", help="Print proposed sublineages to a table.",default=None)
    parser.add_argument("-l", "--labels", help="Print lineage and sample associations to a table formatted for matUtils annotate -c.",default=None)
    parser.add_argument("-t", "--distinction", help="Require that lineage proposals have at least t mutations distinguishing them from the parent lineage or root.",type=int,default=1)
    parser.add_argument("-m", "--minsamples", help="Require that each lineage proposal describe at least m samples.", type=int, default=10)
    parser.add_argument("-w", "--mutweights", help="Path to an optional two (or three) column space-delimited containing mutations and weights (and nodes) to use to weight lineage choices.",default=None)
    parser.add_argument("-g", "--gene", help='Consider only mutations in the indicated gene. Requires that --gtf and --reference be set.', default=None)
    parser.add_argument("-s", "--missense", action='store_true', help="Consider only missense mutations. Requires that --gtf and --reference be set.")
    parser.add_argument("-u", "--cutoff", help="Stop adding serial lineages when at least this proportion of samples are covered. Default 0.95",type=float,default=0.95)
    parser.add_argument("-f", "--floor", help="Minimum score value to report a lineage. Default 0", type=float,default=0)
    parser.add_argument("--gtf", help="Path to a gtf file to apply translation. Use with --reference.")
    parser.add_argument("--reference", help='Path to a reference fasta file to apply translation. Use with --gtf.')
    parser.add_argument("-v","--verbose",help='Print status updates.',action='store_true')
    args = parser.parse_args()
    return args

def main():
    args = argparser()
    t = bte.MATree(args.input)
    mutweights = {}
    if args.gtf != None and args.reference != None:
        if args.verbose:
            print("Performing tree translation and setting weights for mutations based on amino acid changes.")
        translation = t.translate(fasta_file=args.reference,gtf_file=args.gtf)
        for nid, aav in translation.items():
            for aa in aav:
                if args.missense and aa.is_synonymous():
                    continue
                if args.gene != None:
                    if aa.gene != args.gene:
                        continue
                mutweights[(int(aa.nt_index),aa.alternative_nt,nid)] = 1
    if args.mutweights != None:
        mutweights = parse_mutweights(args.mutweights)
    if args.verbose:
        print("Considering {} mutations to have weight.".format(len(mutweights)))
    if args.dump != None:
        dumpf = open(args.dump,'w+')
    if args.clear:
        t.apply_annotations({node.id:[] for node in t.depth_first_expansion()})
    annotes = t.dump_annotations()
    if args.clear:
        assert len(annotes) == 0
    original_annotations = set(annotes.keys())
    global_used_nodes = get_skipset(t, annotes)
    if len(annotes) == 0:
        if args.verbose and not args.clear:
            print("No lineages found in tree; starting from root.")
        annotes = {'L':t.root.id}
    else:
        if args.verbose:
            print("{} annotations found in the tree; identifying candidates for subdivision.".format(len(annotes)))
        annotes = get_outer_annotes(t, annotes)
        if args.verbose:
            print("{} outer annotations found in the tree; identifying sublineages.".format(len(annotes)))
    if args.verbose:
        print("Tree contains {} annotated lineages initially ({} nodes disregarded to prevent retroactive parent assignment).".format(len(annotes),len(global_used_nodes)),file=sys.stderr)
    #keep going until the length of the annotation dictionary doesn't change.
    if args.dump != None:
        print("parent\tparent_nid\tproposed_sublineage\tproposed_sublineage_nid\tproposed_sublineage_score\tproposed_sublineage_size",file=dumpf)
    outer_annotes = annotes
    level = 1
    while True:
        if args.verbose:
            print("Level: ",level)
        new_annotes = {}
        used_nodes = global_used_nodes
        for ann,nid in outer_annotes.items():
            serial = 0
            current_child_lineages = {k:v for k,v in annotes if k.split(".")[:-1] == ann.split(".")}
            print("DEBUG: Found {} child lineages preexisting for lineage {}".format(len(current_child_lineages), ann))
            labeled = set()
            for lin, nid in current_child_lineages.items():
                for s in t.get_leaves_ids(nid):
                    labeled.add(s)
            rbfs = t.breadth_first_expansion(nid, True) #takes the name
            print("DEBUG: Checking annotation {} with {} descendent nodes.".format(nid, len(rbfs)))
            print("DEBUG: Currently, annotation has {} sublineages with {} total samples labeled.".format(len(current_child_lineages),len(labeled)))
            dist_root = dists_to_root(t, t.get_node(nid), mutweights) #needs the node object, not just the name
            while True:
                scdict, leaf_count = get_sum_and_count(rbfs, ignore = labeled, mutweights = mutweights)
                # print("DEBUG: total distances to root {}, total sums {}".format(sum(dist_root.values()),sum([v[0] for v in scdict.values()])))
                best_score, best_node = evaluate_lineage(t, dist_root, nid, rbfs, scdict, args.minsamples, args.distinction, used_nodes)
                if best_score <= args.floor:
                    print("DEBUG: Best doesn't pass threshold with score {} out of {}".format(best_score, args.floor))
                    break
                newname = ann + "." + str(serial)
                while newname in original_annotations:
                    if args.verbose:
                        print("Name {} already in annotation; incrementing".format(newname))
                    serial += 1
                    newname = ann + '.' + str(serial)
                for anc in t.rsearch(best_node.id,True):
                    used_nodes.add(anc.id)
                new_annotes[newname] = best_node.id
                leaves = t.get_leaves_ids(best_node.id)
                if args.dump != None:
                    print("{}\t{}\t{}\t{}\t{}\t{}".format(ann,nid,newname,best_node.id,str(best_score),len(leaves)),file=dumpf)
                for l in leaves:
                    labeled.add(l)
                if len(labeled) >= leaf_count * args.cutoff:
                    break
                serial += 1
                if args.verbose:
                    print("Annotated lineage {} as descendent of {} from level {} with {} descendents".format(newname, ann, level, len(leaves)))
        if not args.recursive:
            annotes.update(new_annotes)
            break
        elif len(new_annotes) == 0:
            break
        else:
            annotes.update(new_annotes)
            outer_annotes = new_annotes
            level += 1
    if args.verbose:
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
        t.apply_annotations(annd)
        t.save_pb(args.output)
    if args.dump != None:
        dumpf.close()
    if args.labels != None:
        labels = {}
        for ann, nid in annotes.items():
            for leaf in t.get_leaves_ids(nid):
                if leaf not in labels:
                    labels[leaf] = [ann]
                else:
                    labels[leaf].append(ann)
        #format this in a way that's parsed by matUtils annotate -c
        with open(args.labels,'w+') as f:
            for l,v in labels.items():
                for ann in v:
                    print("{}\t{}".format(ann,l),file=f)

if __name__ == "__main__":
    main()