'''
Script for the removal of nextclade or other extraneous annotations to prepare a tree for automated Pango lineage annotation.
Nextclade annotations are stored in the first index of standard public trees, so this script overwrites the first index with the second,
except where the second is an automatically maintained lineage, as curated lineages take precedent.
The output tree will contain a single column of annotation containing pango and auto lineages and a blank column.
'''

import bte
import sys
t = bte.MATree(sys.argv[1])
nannd = {}
for n in t.depth_first_expansion():
    if len(n.annotations) >= 2:
        # printer = any(["(" in a for a in n.annotations])
        # if printer:
            # print("ORIGINAL:",n.annotations)
        #overwrite the first column with the value of the second column, then drop the second column
        #unless the second column is an auto annotation AND there's something in the first column, in which case we retain the first column (and still drop the second)
        #also, retain the first column if the second is blank. No need to throw away information. 
        newann = [""]
        #if the second entry is blank, use the first. It may be blank also, which is fine.
        if n.annotations[0] == n.annotations[1]:
            #if you have the same one twice... just collapse it...
            newann[0] = n.annotations[0]
        elif len(n.annotations[1]) == 0:
            #special case- only retain annotations not starting with 1/2 (20A, 20B etc etc) to handle nextclade annotations with no corresponding pango annotation.
            if len(n.annotations[0]) > 0:
                if n.annotations[0][0] not in '12':
                    newann[0] = n.annotations[0]                
        #if the second entry is auto AND the first is not blank, retain the first.
        #also, note that this auto lineage was overridden by a designation if this occurs.
        elif "auto." in n.annotations[1] and len(n.annotations[0]) > 0:
            newann[0] = n.annotations[0]
            print(n.annotations[1], n.annotations[0])
        #otherwise, overwrite the first with it.
        else:
            newann[0] = n.annotations[1]
        #the output only has one entry.
        nannd[n.id] = newann
    elif len(n.annotations) == 1:
        nannd[n.id] = [n.annotations[0]]
    else:
        nannd[n.id] = []
t.apply_annotations(nannd)
t.save_pb(sys.argv[2])