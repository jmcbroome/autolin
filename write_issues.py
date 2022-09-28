import argparse
import pandas as pd
import bte
import numpy as np
import os
from github import Github

def argparser():
    parser = argparse.ArgumentParser(description="Write markdown issues representing the top N proposed sublineages of a report.")
    parser.add_argument("-i", "--input", required=True, help='Path to the lineage report (generate_lineage_report.py) to generate markdown issues from.')
    parser.add_argument("-t", "--tree", required=True, help='Path to the tree-containing protobuf the lineage report was generated from.')
    parser.add_argument("-m", "--metadata", required=True, help='Path to the metadata file corresponding to the tree the report was generate from.')
    parser.add_argument('-n', "--number", default=3, type=int, help="Write markdown reports for the top n lineages.")
    parser.add_argument("-s", "--sort", default='proposed_sublineage_score', help="Choose a column to sort by.")
    parser.add_argument('-p', "--prefix", default='proposal_', help="String to use as prefixes for output files.")
    parser.add_argument("-j", "--jsonsize", default=4000, type=int, help="Maximum size of the json output for each sublineage.")
    parser.add_argument("-l", "--local", action='store_true', help="Set to write issues to local files only. Default posts issues to https://github.com/jmcbroome/auto-pango-designation/issues")
    args = parser.parse_args()
    return args

def write_report(row, prefix):
    fstr = []
    fstr.append("{} is a proposed child lineage of {} including {} samples.".format(row.proposed_sublineage, row.parent, row.proposed_sublineage_size))
    if row.earliest_child != np.nan and row.latest_child != np.nan:
        fstr.append("The earliest sample was found on {} and the latest on {}.".format(row.earliest_child, row.latest_child))
    else:
        fstr.append("Dates could not be identified for these samples.")
    if row.child_regions != np.nan:
        ccount = row.child_regions.count(",") + 1
        common = row.child_regions.split(",")[0]
        commonprop = round(float(row.child_region_percents.split(",")[0])*100,2)
        if ccount == 1:
            fstr.append("It is found in {} only.".format(row.child_regions))
        else:
            fstr.append("It is found in {} countries, most commonly {} where {} of its samples were sequenced.".format(ccount, common, commonprop))
    else:
        fstr.append("Countries could not be identified for these samples.")
    fstr.append("{} has a Bloom Lab escape score of {}, -{} over the parent lineage.".format(row.proposed_sublineage, row.sublineage_escape, row.net_escape_gain))
    spikes = []
    total = 0
    for n in row.aa_changes.split(">"):
        for m in n.split(","):
            if m[0] == 'S':
                spikes.append(m)
            total += 1
    if len(spikes) == 0:
        fstr.append("It is not defined by any spike protein changes. It has {} defining protein changes overall.".format(total))
    else:
        fstr.append("It is defined by the following spike protein changes: {}. There are {} defining protein changes overall.".format(",".join(spikes), total))
    if row.host_jump:
        fstr.append("It represents a zoonotic event!")
    fstr.append("View it on [cov-spectrum]({})".format(row.link))
    return fstr

def write_sample_list(t, mdf, nid, name, prefix):
    with open(prefix + name + "_samples.txt","w+") as outf:
        for s in t.get_leaves_ids(nid):
            try:
                row = mdf.loc[s]
                country = row.country
                date = row.date
            except:
                country = np.nan
                date = np.nan
            print('\t'.join([str(v) for v in [s, country, date]]), file = outf)

def write_json(t, nid, parent_nid, name, prefix, size, metafile = None):
    outn = prefix + name + ".json"
    samples_to_use = t.get_leaves_ids(nid)
    parent_samples = t.get_leaves_ids(parent_nid)
    if len(samples_to_use) + len(parent_samples) > size:
        target = size - len(samples_to_use)
        parents_to_use = np.random.choice(parent_samples, size = target, replace = False)
    else:
        parents_to_use = parent_samples
    if metafile != None:
        t.write_json(outn, samples = samples_to_use + parents_to_use, title = name, metafiles = [metafile])
    else:
        t.write_json(outn, samples = samples_to_use + parents_to_use, title = name, metafiles = [])

def main():
    args = argparser()
    t = bte.MATree(args.tree)
    tn = args.tree.split(".")[0]
    df = pd.read_csv(args.input,sep='\t')
    mdf = pd.read_csv(args.metadata,sep='\t').set_index('strain')
    i = 0
    for ind,d in df.sort_values(args.sort).iterrows():
        print("Recording output for {}".format(d.proposed_sublineage))
        if i >= args.number:
            break
        print("Writing report...")
        report = write_report(d, args.prefix)
        if args.local:
            with open(prefix + d.proposed_sublineage + ".md","w+") as outf:
                print("\n".join(report),file=outf)
        else:
            g = github(os.getenv("API_KEY"))
            r = g.get_user().get_repo("auto-pango-designation")
            titlestring = "Sublineage {} of {}".format(row.proposed_sublineage, row.parent)
            r.create_issue(title=titlestring,body=report)
            
        print("Writing samples...")
        write_sample_list(t, mdf, d.proposed_sublineage_nid, d.proposed_sublineage, args.prefix)
        # submeta_name = d.proposed_sublineage + "_metadata.tsv"
        # print("Writing json...")
        # sdf = mdf.loc[t.get_leaves_ids(d.proposed_sublineage_nid)]
        # print("Submeta extracted; proceeding to write.")
        # if sdf.shape[0] > 0:
            # sdf.reset_index().to_csv(submeta_name,sep='\t',index=False)
        # write_json(t, d.proposed_sublineage_nid, d.parent_nid, d.proposed_sublineage, args.prefix, args.jsonsize, args.metadata)
        # else:
            # write_json(t, d.proposed_sublineage_nid, d.parent_nid, d.proposed_sublineage, args.prefix, args.jsonsize)
        i += 1
    with open(tn+".issues.log","w+") as outf:
        print("Produced files for {} lineage proposals.".format(args.number),file=outf)

if __name__ == "__main__":
    main()