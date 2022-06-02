import bte
import os
import sys
import argparse

#######################################################################################################################################
# This script contains a very early prototype implementation of a single automated lineage labeling algorithm.
# It is not intended to be used for anything other than testing. It works by computing the sum genotype representation 
# of each internal node on the tree and identifying the maximimally representative node, then repeating this process serially
# and hierarchically to generate additional lineages and sublineages. 
#######################################################################################################################################
# MAJOR POINTS OF IMPROVEMENT
# This is a naive implementation and an implementation that uses dynamic programming to track weight contributions in a more 
# informed fashion would likely be substantially faster. 
# Additionally, having a regularization term or informed stopping condition beyond a fixed number of lineages to generate 
# would substantially improve our ability to compare alternative parameters or lineage systems, 
# especially in an hierarchical and online context.
#######################################################################################################################################

def process_mstr(mstr):
    """Read a mutation string and return the chromosome, location, reference, and alternate alleles.
    """
    if ":" in mstr:
        chro = mstr.split(":")[0]
        data = mstr.split(":")[1]
    else:
        chro = None
        data = mstr
    loc = int(data[1:-1])
    ref = data[0]
    alt = data[-1]
    return chro, loc, ref, alt

def compute_representation_scores(nid, tree, halt = None, ranges = []):
    """The representation score of each ancestor of a node is simply the proportion of all mutations
    belonging to a sample that have been accumulated by that ancestor. E.g. if I told you that sample X
    was descended from node N, what percentage of its genotype would you now know?
    """
    mutvec = []
    labels = []
    total = 0
    #To compute the representation score, traverse to the root and record the number of mutations associated with that particular ancestor.
    #ignore mutations which do not fall into the indicated range (if any) and break early if the halt node is reached (so this node does not contribute to nodes above the halt point).
    for n in tree.rsearch(nid,True):
        if n.id == halt:
            break
        mc = 0
        for m in n.mutations:
            if len(ranges) > 0:
                chro, loc, _, _ = process_mstr(m)
                for rchro, rstart, rstop in ranges:
                    if chro == None or rchro == chro:
                        if rstart <= loc <= rstop:
                            mc += 1
                            break
            else:
                mc += 1
        mutvec.append(mc)
        labels.append(n.id)
        total += mc
    if total == 0:
        return [(nid, 0)]
    #once counts have been collected, reverse the count vector, compute a proportion accumulation vector, and return the proportion accumulation vector zipped with the reversed label vector.
    props = []
    accumulated = 0
    for mv in mutvec[::-1]:
        accumulated += mv
        props.append(accumulated / total)
        
    return zip(labels[::-1],props)

def compute_representation_map(tree, skip = set(), subroot = None, ranges = [], sweights = {}):
    """Compute the representation map for all samples on the tree. Scores are computed per sample.
    Scores are summed for each internal node from each sample contributing to that node. 
    """
    repmap = {}
    #if computing for a subtree, start at the indicated node and halt weight computation above the indicated node.
    if subroot == None:
        base = tree.root.id
    else:
        base = subroot
    for l in tree.get_leaves_ids(base):
        #when computing lineages serially, skip samples that were already assigned to another lineage.
        if l not in skip:
            #if this sample has a special weight multiplier, apply it. 
            mult = sweights.get(l,1)
            for anc, prop in compute_representation_scores(l, tree, subroot, ranges):
                if anc not in repmap:
                    repmap[anc] = 0
                repmap[anc] += prop * mult
    return repmap

def label_lineages_serial(tree, count = 5, subroot = None, coverage_threshold = 0.95, ranges = [], sweights={}):
    """To compute multiple lineage labels serially (meaning one is not a subset of another and all labels are mutually exclusive)
    We repeatedly compiute and assign the maximally representative node while ignoring samples which are descended from a previously assigned node.
    This can be repeated until a fixed number of lineages is reached or a total percentage of samples are covered by at least one annotation.
    """
    labels = {}
    samples_labeled = set()
    if subroot == None:
        subroot = tree.root.id
    total_sample_count = len(tree.get_leaves_ids(subroot))
    labels_generated = 0
    while len(samples_labeled) < (total_sample_count * coverage_threshold) and labels_generated < count:
        repmap = compute_representation_map(tree, samples_labeled, subroot, ranges, sweights)
        if len(repmap) > 0:
            target = max(repmap.items(),key=lambda x:x[1])[0]
            labels[target] = str(labels_generated)
                
            for l in tree.get_leaves_ids(target):
                samples_labeled.add(l)
        labels_generated += 1
    
    return labels

def label_lineages_hierarchical(tree, serial_count = 3, levels = 2, coverage_threshold = 0.95, ranges=[], size=0, sweights={}, current_annotations = {}):
    """To compute multiple lineages labels hierarchically (meaning one label is a subset or child of another) we compute a set of labels
    serially, then break the tree into subtrees at each of those serial labels and compute an additional set of serial labels independently
    from each subtree. This is repeated until a fixed number of levels have been computed.
    """
    #hierarchies are defined as a dictionary of dictionaries of labels
    level_labels = {}
    if len(current_annotations) != 0:
        level_labels[-1] = current_annotations
    for level in range(levels):
        if len(level_labels) == 0:
            labels = label_lineages_serial(tree, count=serial_count, coverage_threshold=coverage_threshold, ranges=ranges, sweights=sweights)
            level_labels[level] = labels
        else:
            if level not in level_labels:
                level_labels[level] = {}
            for nid,outer_label in level_labels[level-1].items():
                leafcount = len(tree.get_leaves_ids(nid))
                if leafcount >= size:
                    labels = label_lineages_serial(tree, count=serial_count, subroot=nid, coverage_threshold=coverage_threshold, ranges=ranges, sweights=sweights)
                    level_labels[level].update({k:outer_label+"."+v for k,v in labels.items()})
    return level_labels

def parse_current_annotations(af):
    """Parse a set of annotations to start with from a file (tab delimited, node lineage) and return a dictionary of the form {nid:label}
    """
    annotations = {}
    with open(af) as inf:
        for entry in inf:
            nid, label = entry.strip().split("\t")
            if nid != "node_id":
                annotations[nid] = label
            else:
                print("WARNING: Skipping current annotation entry 'node_id'")
    return annotations

def parse_ranges(bf):
    """Parse a list of ranges from a bed format file (tab delimited, chromosome start stop). Ranges are inclusive.
    """
    ranges = []
    with open(bf) as inf:
        for entry in inf:
            chro, start, stop = entry.strip().split()
            ranges.append((chro,int(start),int(stop)))
    return ranges
    
def parse_weights(wf):
    """Parse a set of sample weights to use from a file (tab delimited, sample weight). Samples not included in the file are assigned weight 1.
    """
    weightd = {}
    with open(wf) as inf:
        for entry in inf:
            sample, weight = entry.strip().split("\t")
            weightd[sample] = float(weight)
    return weightd

def argparser():
    parser = argparse.ArgumentParser(description="Generate lineage assignments for a phylogeny based on total genotype representation.")
    parser.add_argument("-i", "--input", required=True, help='Path to protobuf to annotate.')
    parser.add_argument("-a", "--count", type=int, help="Maximum number of annotations to use for each annotation level. Default is 10.", default = 10)
    parser.add_argument("-l", "--levels", type=int, help="Maximum number of levels of hierarchical annotation to perform. Default 1.", default = 1)
    parser.add_argument("-z", "--size", type=int, help="Maximum size for a final annotation level. Lineages larger than this will be further subdivided as maximum levels allow. Default behavior does not use this value.",default=0)
    parser.add_argument("-o", "--output", help="Name of the file to write lineage node assignments to.", default = None)
    parser.add_argument("-m", "--metadata", help="Name of the file to write samples and their highest-level annotation to, for use with Nextstrain.", default = None)
    parser.add_argument("-t", "--threshold", type=float, help="Set to a percentage of all samples to cover with annotations in each level. Annotations will generate until either this threshold is reached or the indicated number with -a is generated. Default 0.95", default=0.95)
    parser.add_argument("-s", "--silent", action='store_true', help="Supress stderr output from this tool. Note that this may silence error messages.")
    parser.add_argument("-r", "--range", help="Pass a bed file containing intervals; mutations within these intervals will be used to construct genotypes. If this option is unused, the whole genome is used.", default=None)
    parser.add_argument("-w", "--weights", help="Pass a two-column tab-delimited file containing sample names in the first column and weights in the second column. Default for any sample not included is weight 1.", default=None)
    parser.add_argument("-c", "--current", help="Pass a tab-delimited file containing node ids and their corresponding annotations. These will be used to start the first level of annotations.", default=None)
    args = parser.parse_args()
    return args

def main():
    args = argparser()
    if args.silent:
        f = open(os.devnull, 'w')
        sys.stderr = f
    ranges = []
    if args.range != None:
        ranges = parse_ranges(args.range)
    ad = {}
    if args.current != None:
        ad = parse_current_annotations(args.current)
    wd = {}
    if args.weights != None:
        wd = parse_weights(args.weights)
    t = bte.MATree(args.input)
    labeld = label_lineages_hierarchical(t, serial_count=args.count,levels=args.levels, coverage_threshold=args.threshold,ranges=ranges,size=args.size,sweights=wd,current_annotations=ad)
    if args.output != None:
        with open(args.output,"w+") as outf:
            print("node_id\tlabel",file=outf)
            for l, ld in labeld.items():
                for nid, nidlabel in ld.items():
                    print(nid+"\t"+nidlabel,file=outf)
    if args.metadata != None:
        with open(args.metadata,"w+") as outf:
            labeled = set()
            print("sample\tlineage",file=outf)
            for level in reversed(range(args.levels)):
                for k,v in labeld[level].items():
                    for l in t.get_leaves_ids(k):
                        if l not in labeled:
                            labeled.add(l)
                            print(l+"\tlin_"+str(v),file=outf)

if __name__ == "__main__":
    main()