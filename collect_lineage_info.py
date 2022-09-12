import bte
# import numpy as np
import datetime as dt
from os.path import exists
import sys

cdate = sys.argv[1]

def get_desc_date(s):
    try:
        datestr = dt.datetime.strptime(s.split("|")[-1],"%Y-%m-%d")
        return datestr
    except:
        return None

def get_target_files(date):
    if type(date) == str:
        datestr = date
    else:
        datestr = date.strftime("%Y-%m-%d")
    return (datestr + ".pango.txt", "run_1/" + datestr + ".proposed.txt", "daily_trees/public-" + datestr + ".all.masked.nextclade.pangolin.pb.gz")

def get_new_annotations(annf, seen = set()):
    new = {}
    with open(annf) as inf:
        for entry in inf:
            spent = entry.strip().split('\t')
            if spent[0] == 'clade':
                continue
            if spent[0] not in seen:
#                 print(annf,entry)
                new[spent[0]] = spent[1]
                seen.add(spent[0])
    return new, seen

pangof,propf,treef = get_target_files(cdate)
seen = set()
if exists('tmp_seen.txt'):
    with open("tmp_seen.txt") as inf:
        for c in inf:
            seen.add(c.strip())
    assert len(seen) > 0

if exists(pangof):
#         pango = pd.read_csv(pangof,sep='\t')
#         pango.set_index("clade",inplace=True)
    new_anns,seen = get_new_annotations(pangof,seen)
    if len(new_anns) == 0:
        print("No new annotations on {}; continuing".format(cdate.strftime("%Y-%m-%d")))
        # continue
        exit(0)
    assert len(seen) > 0
    with open("tmp_seen.txt","w+") as outf:
        for c in seen:
            print(c,file=outf)
else:
    print("Pango file for {} missing; continuing".format(cdate.strftime("%Y-%m-%d")))
    # continue
    exit(0)
#proposed = pd.read_csv(propf,sep='\t')
t = bte.MATree(treef)
# outd = {k:[] for k in ['Lineage','Parent','DateAdded','EarliestDate','LatestDate','Count']}
outf = open(cdate+".linfo.tsv","w+") 
print('\t'.join(['DateAdded','Lineage','Parent','EarliestDate','LatestDate','Count']),file=outf)
for c, rid in new_anns.items():
#         print(c,rid)
    node = t.get_node(rid)
    try:
        parent_lineages = ",".join(node.parent.most_recent_annotation())
    except AttributeError:
        parent_lineages = 'None'
    descendents = t.get_leaves_ids(rid)
#         track = date_tracker.get(c,None)
#         if track == None:
    descendent_dates = [get_desc_date(s) for s in descendents if get_desc_date(s) != None]
    if len(descendent_dates) == 0:
        print("Could not identify dates for samples for node {} on day {}!".format(rid,cdate))
        edate = None
        mdate = None
    else:
        edate = min(descendent_dates).strftime("%Y-%m-%d")
        mdate = max(descendent_dates).strftime("%Y-%m-%d")
    print("\t".join([str(v) for v in [cdate, c, parent_lineages, edate, mdate, len(descendents)]]),file=outf)
outf.close()
