import bte
import sys
import argparse
import yaml

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

def process_metadata(metadata):
    """Read metadata file and returns a dict using the sample (leaf) as the key and the other fields(date, region, etc) in a list as the value.
    """
    ans = []
    metadata_dict = {}
    # open .tsv file
    with open(metadata) as f:
        # Read data line by line
        for line in f:
            # split data by tab
            # store it in list
            l=line.split('\t')
            # append list to ans
            ans.append(l)
            metadata_dict[l[0]] = l[1:-1]
    return metadata_dict

def pango_lineage_naming(annote, ser, alias_key, desc_depth):
    #### desc_depth (4 per pango guidelines) can be altered by user. This integer specifies how many levels of descendants constitues a letter change for the lineage name
    #Throw Away Letters = "I", "O", "X" per pango guidelines
    throw_away_letters = []
    throw_away_letters.append(ord("X"))   # for recombinant lineages, ignored for now
    throw_away_letters.append(ord("I"))
    throw_away_letters.append(ord("O"))

    proposed_name = annote + "." + str(ser)
    descendant_count = proposed_name.count('.')
    alias = proposed_name
    proposed_name = list(proposed_name)
    count  = 0
    while (proposed_name[count] != '.'):
        if proposed_name[count].isdigit():    #ex 21D (eta).x.x.x   essentially ignoring nextstrain info also not including in alias key
            proposed_name = ''.join(proposed_name)
            return str(proposed_name), alias_key
        letter = ord(proposed_name[count])
        if letter in throw_away_letters:
            proposed_name[count] = chr(letter + 1)
        count += 1
    if descendant_count >= desc_depth:
        char = ''
        count = 0
        beginning_letters = []
        while (proposed_name[count] != '.'):
            char = proposed_name[count]
            beginning_letters.append(char)
            count += 1
        k = -1
        for i in range(1, len(beginning_letters)):
            if (beginning_letters[k] == 'Z') or (beginning_letters[k] == 'z'):
                if (len(beginning_letters) > k*(-1)):   #if current letter is at least not the last letter
                    # CZB.3.2.2.1 -> CZC.1        CZZ.3.2.2.1 -> DAA.1    CBD.3.2.2.2 -> CEB.2  CBZ.3.2.2.2 -> CCA.2
                    beginning_letters[k] = "A"
                    if (beginning_letters[k-1] != "Z") and (beginning_letters[k-1] != "z"):
                        letter = ord(beginning_letters[k - 1])  #increments letter before if not z
                        beginning_letters[k - 1] = chr(letter + 1)
                    k -= 1
            else:   # increment normal
                letter = ord(beginning_letters[k])
                beginning_letters[k] = chr(letter + 1)
                # beginning_letters is incremented and ready now.
                break
        dot_count = 0
        charcount = 0
        for char in proposed_name:
            charcount += 1
            if char == '.':
                dot_count += 1
            if dot_count >= 4:
                break
        proposed_name = ''.join(beginning_letters + proposed_name[charcount-1:])
        len(beginning_letters)
    proposed_name = ''.join(proposed_name)
    alias_key[proposed_name] = alias
    return str(proposed_name), alias_key

def generate_yaml_report(nameid, node, alias_key, md_dict, tree):   #call for every lineage name run (whenever pango_lineage_naming is called)
    """
    ####EXAMPLE####
    name: B.1
    unaliased_name: A.1.1.1
    parent: A.1.1
    designation_date: "2022-06-10"
    defining_snps:
        - pos: 77383
        nucleotide: A
    reference_sequences:
        - source: genbank
        accession: ON563414
        isolate: MPXV_USA_2022_MA001"""
    # not working as of now. The concept is to grab the parsed information from the metadata file (md_dict) which uses the sample name (leaf)
    # as the key and all other fields (date, region, etc) are in a list as the value. The idea would be to write the needed information from the
    # md_dict into a yaml file for the specific lineage being looked at currently. It should be noted that figuring out the earliest and latest
    # designation dates for the lineage requires looking at all leaves of a lineage and keeping track of earliest and latest leaves (along with
    # the nodes associated with them). All other information such as accession, isolate, pos (position) is derived from the latest sample but
    # this can be changed, it made sense to be the latest sample's information though.

    # include date of earliest lineage detection sample, latest lineage detection sample, and current date of protobuf
    
    if nameid in alias_key:
        curr_alias = alias_key[nameid]
    else:
        return
    # need to assess each child sample (leaf) to pick the correct info to sum up the majority of leaves. all from parsed metadata file
    child_list = []
    for children in node.children:    # gets all the leaves for the specific lineage being looked at
        child_list.append(children.id)

    earliest_date = 0   # cannot be 0
    latest_date = 0
    earliest_date_sample = ""
    latest_date_sample = ""
    refined_child_list = []
    for leaf in child_list:
        if leaf in md_dict:
            refined_child_list.append(leaf)
    for leaf in refined_child_list:
        if leaf in md_dict:
            date = md_dict[leaf][1]
            if date.count("-") == 0:
                if int(date) > latest_date:
                    # print("CHANGED 1.0")
                    latest_date = int(date)
                    latest_date_sample = leaf
                elif int(date) < earliest_date:
                    # print("CHANGED 1.1")
                    earliest_date = int(date)
                    earliest_date_sample = leaf
            else:
                date = date.split("-")[0]
                if int(date) > latest_date:
                    # print("CHANGED 2.0")
                    latest_date = int(date)
                    latest_date_sample = leaf
                elif int(date) < earliest_date:
                    # print("CHANGED 2.1")
                    earliest_date = int(date)
                    earliest_date_sample = leaf
    parent = nameid.split('.')[:-1]
    # print("PARENT, NODE: ")
    # print(parent, nameid)
    # print("earliest date sample: ")
    # print(earliest_date_sample)
    # print(md_dict[earliest_date_sample])
    report_list = {
    # "Shopping List": OrderedDict({
    "name": nameid,
    "unaliased_name": curr_alias,
    "parent": parent,
    "earliest_designation_date": md_dict[earliest_date_sample][1],
    "latest_designation_date": md_dict[latest_date_sample][1],
    "defining_snps": {
        "- pos": md_dict[latest_date_sample][5],
        "  nucleotide": tree.get_haplotype(latest_date_sample),
        },
    "reference_sequences": {
        "- source": "genbank",
        "  accession": md_dict[latest_date_sample][0],
        "  isolate": md_dict[latest_date_sample][2]
        }
    }

    with open(str(node.id) + '.yml', 'w') as yaml_file:
        yaml.dump(report_list, yaml_file, default_flow_style=False)
    return

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
    parser.add_argument("-md", "--metadata", help='Path to a metadata tsv file.')
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
    if args.metadata != None:
        metadata = process_metadata(args.metadata)
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
    alias_key = {}
    while True:
        if args.verbose:
            print("Level: ",level)
        new_annotes = {}
        used_nodes = global_used_nodes
        for ann,nid in outer_annotes.items():
            serial = 0
            current_child_lineages = {k:v for k,v in annotes.items() if k.split(".")[:-1] == ann.split(".") and k != ann}
            labeled = set()
            for lin, cnid in current_child_lineages.items():
                for s in t.get_leaves_ids(cnid):
                    labeled.add(s)
            # if len(current_child_lineages) > 0:
                # print("DEBUG: Found {} child lineages preexisting for lineage {}; {} samples prelabeled from {} total samples".format(len(current_child_lineages), ann, len(labeled), len(t.get_leaves_ids(nid))))
            rbfs = t.breadth_first_expansion(nid, True) #takes the name
            # print("DEBUG: Checking annotation {} with {} descendent nodes.".format(nid, len(rbfs)))
            dist_root = dists_to_root(t, t.get_node(nid), mutweights) #needs the node object, not just the name
            while True:
                scdict, leaf_count = get_sum_and_count(rbfs, ignore = labeled, mutweights = mutweights)
                # print("DEBUG: total distances to root {}, total sums {}".format(sum(dist_root.values()),sum([v[0] for v in scdict.values()])))
                best_score, best_node = evaluate_lineage(t, dist_root, nid, rbfs, scdict, args.minsamples, args.distinction, used_nodes)
                if best_score <= args.floor:
                    # print("DEBUG: Best doesn't pass threshold with score {} out of {}".format(best_score, args.floor))
                    break
                newname = ann + "." + str(serial)

                while newname in original_annotations:
                    serial += 1
                    newname = ann + '.' + str(serial)
                for anc in t.rsearch(best_node.id,True):
                    used_nodes.add(anc.id)
                new_annotes[newname] = best_node.id
                leaves = t.get_leaves_ids(best_node.id)
                if args.dump != None:
                    name, alias_key = pango_lineage_naming(ann, serial, alias_key, 4)  # renames lineage if needed per pango guidelines
                    print("{}\t{}\t{}\t{}\t{}\t{}".format(name,nid,name + '.' + str(serial),best_node.id,str(best_score),len(leaves)),file=dumpf)
                    # generate_yaml_report(name, best_node, alias_key, metadata, t)    # for yaml file writing
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
    print("ALIAS KEY: ")
    print(alias_key)

if __name__ == "__main__":
    main()