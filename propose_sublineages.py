from curses.ascii import isdigit
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

def pango_lineage_naming(annote, ser, alias_key):   # may need to pass in a list of all lineage names to compare and create an accurate new name 
    # mostly for adding to existing pango lineage names
    # no hardcoding names. must be modular since it will be expanded for other disease

    # 1 and 15 for new parameters
    # print("ANNOTE: ", annote)

    #QUESTION: should greek letters be converted to pango lineage naming convention? What would they be incremented to? 

    #should change name even if it doesnt have 4 levels of descendants if the name is already incorrect. Check for throw away letters or weird symbols
    # need to figure out the correct way to change these


    ###### should aim to have # of decimal places that cause increment of letter to be passed into the function (its 4 normally)
    ###### this should also mean that the incrementing algorithm should be in a loop that runs with less hard-coded steps, have to think about this more



    # create alias that could be on the dump file that tells you the summary alias. 
    # "BD" = "B.1.1.529.1.17.2"  for example. just so no information is lost and it can be translated back from summary.

    # following pango naming rules except for recombination lineages (X)

    # omitted "I", "O", "X" per pango guidelines (section 2.1.b)

    
    letters = {"0":"A", "1":"B", "2":"C", "3":"D", "4":"E", "5":"F", "6":"G", "7":"H", "8":"J", "9":"K", "10":"L", "11":"M", "12":"N", "13":"P",
    "14":"Q", "15":"R", "16":"S", "17":"T", "18":"U", "19":"V", "20":"W", "21":"Y", "22":"Z"}


    #Throw Away Letters = "I", "O", "X"
    throw_away_letters = []
    throw_away_letters.append(ord("X"))
    throw_away_letters.append(ord("I"))
    throw_away_letters.append(ord("O"))

    # encountering bug in dump file where some names of lineages have numbers in front such as 20G.13.0.0 or 21B (Kappa).2 or XB.0
    # thinking this may be involving how im incrementing letters in ascii and converting back.
    # this could possibly be due to inconsistency in the pango naming algorithm logic. Will change this.

    proposed_name = annote + "." + str(ser)
    descendant_count = proposed_name.count('.')
    alias = proposed_name
    proposed_name = list(proposed_name)
    
    # loop to check if throw away letters are already present in the current name and change them. ( or weird symbols ex: 21D (Eta))

    count  = 0
    while (proposed_name[count] != '.'):
        # currlett = roposed_name[count]
        if proposed_name[count].isdigit():    #ex 21D (eta).x.x.x   essentially ignoring nextstrain info also not including in alias key
            proposed_name = ''.join(proposed_name)
            return str(proposed_name), alias_key
        letter = ord(proposed_name[count])
        if letter in throw_away_letters:
            proposed_name[count] = chr(letter + 1)
        count += 1

    # if (proposed_name[0] == 'L') and (proposed_name[1] == '.'):
    #     numb = proposed_name[2]   #doesnt account for double digit numbers    must fix this

    #     if numb in letters:
    #         letter = letters[numb]
    #         # print("LENG: " + str(len(proposed_name))
    #         proposed_name.pop(0)
    #         proposed_name.pop(1)
    #         proposed_name.pop(2)
    #         proposed_name.append(letter)


    #     # ###### maybe shouldnt be in else
    #     # else:  #number is bigger than 22   number%22 until less than 22? that is remainder that will be converted into letter
    #     #     remain = 23
    #     #     while (remain > 22):
    #     #         remain = numb%22
    #     # #####

    #could be if descendant_count >= inputted descendant depth number
    if descendant_count >= 4:  #too many periods for descendants therefore should increment the alphabetical name up 1 letter

        # Question:  what does A.1.3.4.28 turn into? B.28? this would lose the 1.3.4 information since its not 1.1.1  = bc it is the same however alias key will give more info

    #     #if [0] is Z and [1] is '.', make AA
    #     # if [0] is A and [1] is Z, make BA and so on
    #     # if [1] is not '.' (meaning its a letter), increment second letter (unless Z then do AA)
    #     # if [0] is Z and [1] is Z, next is AAA
        if (proposed_name[1] == '.'):
            # regular increment to next letter with descendant number after  ex: A.1.1.1.28 -> B.28
            if ((proposed_name[0] != 'Z') and (proposed_name[0] != 'z')):
                # print("Proposed name0: ", proposed_name)
                letter = ord(proposed_name[0])
                if letter in throw_away_letters:
                    proposed_name[0] = chr(letter + 2)
                else:
                    proposed_name[0] = chr(letter + 1)
                dot_count = 0
                count = 0
                # desc_numb = 0
                for char in proposed_name:
                    count += 1
                    if char == '.':
                        dot_count += 1
                    if dot_count >= 4:
                        # desc_numb = proposed_name[count+1]
                        break

                new_prop_name = []
                # print("Proposed name1: ", proposed_name)
                for i in range(0, 1):  # wouldn't account for numbers after 4th decimal place. ####Must fix this.#####
                    # print("I: ", i, " proposed_name[i] = ", proposed_name[i])
                    new_prop_name.append(proposed_name[i])    #error on full data set
                # print("Proposed name2: ", proposed_name)
                proposed_name = ''.join(new_prop_name + proposed_name[count-1:])
                # print("Proposed name2: ", proposed_name)
                # for i in range(1, count - 1):
                #     proposed_name.pop(i)
                alias_key[proposed_name] = alias
                
            # if [0] is Z and [1] is '.', make AA  ex: Z.1.1.1.53 -> AA.53
            elif ((proposed_name[0] == 'Z') or (proposed_name[0] == 'z')):  # if first letter is already Z
                dot_count = 0
                count = 0
                # desc_numb = 0
                for char in proposed_name:
                    count += 1
                    if char == '.':
                        dot_count += 1
                    if dot_count >= 4:
                        # desc_numb = proposed_name[count+1]
                        break
                new_prop_name = []        
                # for i in range(0, count - 1):
                #     proposed_name.pop(i)
                new_prop_name.insert(0, "A")   
                new_prop_name.insert(0, "A")  #should now be AA.x   where x is some descendant number 
                proposed_name = ''.join(new_prop_name + proposed_name[count-1:])
                alias_key[proposed_name] = alias

        # if [0] is A and [1] is Z, make BA and so on
        else:  # first period is not in 1st index. dealing with at least a double character lineage name
            if (proposed_name[2] == '.'):
                # if [0] is some letter and [1] is not Z and [2] is '.' (2 letters) increment [1] to next letter. ex: BC.1.1.1.33 -> BD.33
                if ((proposed_name[1] != 'Z') and (proposed_name[1] != 'z')):
                    # print("Proposed name0: ", proposed_name)   # at beginning of 
                    letter = ord(proposed_name[1])
                    if letter in throw_away_letters:
                        proposed_name[1] = chr(letter + 2)
                    else:
                        proposed_name[1] = chr(letter + 1)
                    dot_count = 0
                    count = 0
                    # desc_numb = 0
                    for char in proposed_name:
                        count += 1
                        if char == '.':
                            dot_count += 1
                        if dot_count >= 4:
                            # desc_numb = proposed_name[count+1]
                            break
                    # print("Proposed name: ", proposed_name)
                    new_prop_name = []
                    for i in range(0, 2):
                        # print("I: ", i, " proposed_name[i] = ", proposed_name[i])
                        new_prop_name.append(proposed_name[i])    #error on full data set
                    # print("Proposed name2: ", proposed_name)
                    proposed_name = ''.join(new_prop_name + proposed_name[count-1:])
                    # print("Proposed name3: ", proposed_name)
                    alias_key[proposed_name] = alias
                
                # if [0] is some letter and [1] is Z and [2] is '.' (2 letters) increment [1] to next letter unless already Z.  ex: BZ.1.1.1.51 -> CA.51
                elif ((proposed_name[1] == 'Z') or (proposed_name[1] == 'z')) and ((proposed_name[0] != 'Z') and (proposed_name[0] != 'z')):  # if first letter is already Z
                    letter = ord(proposed_name[0])
                    # print("PROPOSED NAME 1: ", proposed_name)
                    if letter in throw_away_letters:
                        proposed_name[0] = chr(letter + 2)
                    else:
                        proposed_name[0] = chr(letter + 1)
                    
                    dot_count = 0
                    count = 0
                    # desc_numb = 0
                    for char in proposed_name:
                        count += 1
                        if char == '.':
                            dot_count += 1
                        if dot_count >= 4:
                            # desc_numb = proposed_name[count+1]
                            break
                    new_prop_name = []
                    # for i in range(1, count - 1):
                    #     proposed_name.pop(i)
                    new_prop_name.append(proposed_name[0])
                    new_prop_name.append("A")
                    # new_prop_name.insert(1, "A")   #should now be XA.x   where x is some descendant number and X is some letter
                    proposed_name = ''.join(new_prop_name + proposed_name[count-1:])
                    # print("PROPOSED NAME 2: ", proposed_name)
                    alias_key[proposed_name] = alias

                # if [0] is Z and [1] is Z and [2] is '.' (2 letters) change into AAA.  ex: ZZ.1.1.1.51 -> AAA.51
                elif ((proposed_name[1] == 'Z') or (proposed_name[1] == 'z')) and ((proposed_name[0] == 'Z') or (proposed_name[0] == 'z')):  # if first letter is already Z
                    dot_count = 0
                    count = 0
                    # desc_numb = 0
                    for char in proposed_name:
                        count += 1
                        if char == '.':
                            dot_count += 1
                        if dot_count >= 4:
                            # desc_numb = proposed_name[count+1]
                            break
                    new_prop_name = []
                    # for i in range(0, count - 1):
                    #     proposed_name.pop(i)
                    new_prop_name.insert(0, "A")
                    new_prop_name.insert(0, "A")
                    new_prop_name.insert(0, "A")   #should now be XA.x   where x is some descendant number and X is some letter
                    new_prop_name = ''.join(new_prop_name + proposed_name[count-1:])
                    alias_key[proposed_name] = alias

            elif (proposed_name[3] == '.'):
                # if [0] is some letter and [1] is some letter and [2] is not Z and [3] is '.' (3 letters) increment [2].  ex: AAA.2.1.3.24 -> AAB.24
                if (proposed_name[2] != 'Z') and (proposed_name[2] != 'z'):

                    letter = ord(proposed_name[2])
                    if letter in throw_away_letters:
                        proposed_name[2] = chr(letter + 2)
                    else:
                        proposed_name[2] = chr(letter + 1)
                    # proposed_name[0] = chr(letter + 1)
                    dot_count = 0
                    count = 0
                    # desc_numb = 0
                    for char in proposed_name:
                        count += 1
                        if char == '.':
                            dot_count += 1
                        if dot_count >= 4:
                            # desc_numb = proposed_name[count+1]
                            break
                    
                    new_prop_name = []
                    # print("Proposed name1: ", proposed_name)
                    for i in range(0, 3):  # wouldn't account for numbers after 4th decimal place. ####Must fix this.#####
                        # print("I: ", i, " proposed_name[i] = ", proposed_name[i])
                        new_prop_name.append(proposed_name[i])    #error on full data set
                    # print("Proposed name2: ", proposed_name)
                    proposed_name = ''.join(new_prop_name + proposed_name[count-1:])
                    # print("PROP NAME: ", proposed_name)
                    # for i in range(3, count - 1):
                    #     proposed_name.pop(i)
                    alias_key[proposed_name] = alias
                
                # if [0] is not Z and [1] is not Z and [2] is Z and [3] is '.' (3 letters) increment [2].  ex: AAZ.1.2.2.27 -> ABA.27
                elif ((proposed_name[2] == 'Z') or (proposed_name[2] == 'z')) and ((proposed_name[1] != 'Z') and (proposed_name[1] != 'z')):
                    letter = ord(proposed_name[1])
                    if letter in throw_away_letters:
                        proposed_name[1] = chr(letter + 2)
                    else:
                        proposed_name[1] = chr(letter + 1)
                    dot_count = 0
                    count = 0
                    # desc_numb = 0
                    for char in proposed_name:
                        count += 1
                        if char == '.':
                            dot_count += 1
                        if dot_count >= 4:
                            # desc_numb = proposed_name[count+1]
                            break
                    new_prop_name = []
                    # for i in range(2, count - 1):
                    #     proposed_name.pop(i)
                    new_prop_name.append(proposed_name[0])
                    new_prop_name.append(proposed_name[1])
                    new_prop_name.append("A")
                    proposed_name = ''.join(new_prop_name + proposed_name[count-1:])
                    alias_key[proposed_name] = alias

                # if [0] is not Z and [1] is Z and [2] is Z and [3] is '.' (3 letters) increment [2].  ex: AZZ.1.2.2.27 -> BAA.27
                elif ((proposed_name[2] == 'Z') or (proposed_name[2] == 'z')) and ((proposed_name[1] == 'Z') or (proposed_name[1] == 'z')) and ((proposed_name[0] != 'Z') and (proposed_name[0] != 'z')):
                    letter = ord(proposed_name[0])
                    if letter in throw_away_letters:
                        proposed_name[0] = chr(letter + 2)
                    else:
                        proposed_name[0] = chr(letter + 1)
                    dot_count = 0
                    count = 0
                    # desc_numb = 0
                    for char in proposed_name:
                        count += 1
                        if char == '.':
                            dot_count += 1
                        if dot_count >= 4:
                            # desc_numb = proposed_name[count+1]
                            break
                    new_prop_name = []
                    # for i in range(1, count - 1):
                    #     proposed_name.pop(i)
                    new_prop_name.append(proposed_name[0])
                    new_prop_name.append("A")
                    new_prop_name.append("A")
                    proposed_name = ''.join(new_prop_name + proposed_name[count-1:])
                    alias_key[proposed_name] = alias


            # if (proposed_name != 'Z') or (proposed_name != 'z'):

            #     proposed_name[0] = chr(letter + 1)
            #     str(proposed_name)
            # else:
            #     proposed_name.insert(1, 'A')
    
        ### UNCOMMENT START HERE 
        ##### new faster method that would allow descendant depth number to be inputted as any number (not just 4)
        # ####################################################
        # #working on more condensed version of naming algorithm
        # # must figure out where first decimal place is:
        # char = ''
        # count = 0
        # beginning_letters = []
        # while (proposed_name[char] != '.'):
        #     char = proposed_name[count]
        #     beginning_letters.append(char)
        #     count += 1
        # # at this point the count is where the decimal place is

        # #algorithm takes place here:
        
        # #smallest (most specific) character

        # # for i in range(1, len(beginning_letters)):
        
        # # k = -i
        # k = -1

        # # while ((beginning_letters[k] != 'Z') and (beginning_letters[k] != 'z')):   #when to stop? when the current letter is not z anymore.   # wb ZCZ -> ZDA
        # for i in range(1, len(beginning_letters)):   # figure out break condition or maybe this works fine

        #     if (beginning_letters[k] == 'Z') or (beginning_letters[k] == 'z'):
        #         if (len(beginning_letters) > k*(-1)):   #if current letter is at least not the last letter
        #             # if (beginning_letters[-2]
        #             # k = -i -1
        #             k -= 1
        #             # beginning_letters[-i - 1]


        #         else:
                    

        #             continue

        #     else:   # increment normal

        #         letter = ord(beginning_letters[k])

        #         beginning_letters[k] = chr(letter + 1)

        #         # beginning_letters is incremented and ready now.
        #         break


        # # dot_count = 0
        # # count = 0
        # # # desc_numb = 0
        # # for char in proposed_name:
        # #     count += 1
        # #     if char == '.':
        # #         dot_count += 1
        # #     if dot_count >= 4:
        # #         # desc_numb = proposed_name[count+1]
        # #         break

        # # for i in range(1, count - 1):
        # #     proposed_name.pop(i)
        # # alias_key[proposed_name] = alias


    
        # len(beginning_letters)
        
        # ##################################################
        #### UNCOMMENT END HERE

    proposed_name = ''.join(proposed_name)
    # print("name: ", proposed_name)
    alias_key[proposed_name] = alias
    return str(proposed_name), alias_key


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

def evaluate_lineage(t, dist_to_root, anid, candidates, sum_and_count, minimum_size = 0, minimum_distinction = 0):
    """Evaluate every descendent branch of lineage a to propose new sublineages.

    Args:
        t (MATree): The tree.
        a (str): The lineage annotation node to check.
    """
    good_candidates = []
    for c in candidates:
        if not c.is_leaf():
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
    alias_key = {}
    while True:
        print("Level: ",level)
        new_annotes = {}
        for ann,nid in outer_annotes.items():
            serial = 0
            labeled = set()
            rbfs = t.breadth_first_expansion(nid, True) #takes the name
            dist_root = dists_to_root(t, t.get_node(nid), mutweights) #needs the node object, not just the name
            while True:
                scdict, leaf_count = get_sum_and_count(rbfs, ignore = labeled, mutweights = mutweights)
                best_score, best_node = evaluate_lineage(t, dist_root, nid, rbfs, scdict, args.minsamples, args.distinction)
                if best_score <= 0:
                    break
                newname = ann + "." + str(serial)
                new_annotes[newname] = best_node.id
                if args.dump != None:
                    print("ADDING LINEAGES FROM LEVEL {}".format(level))
                    name, alias_key = pango_lineage_naming(ann, serial, alias_key)
                    print("{}\t{}\t{}\t{}\t{}".format(name,nid,name + '.' + str(serial),best_node.id,str(best_score)),file=dumpf)
                    # print("{}\t{}\t{}\t{}\t{}".format(ann,nid,newname,best_node.id,str(best_score)),file=dumpf)

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
    print("ALIAS KEY: ")
    print(alias_key)

if __name__ == "__main__":
    main()