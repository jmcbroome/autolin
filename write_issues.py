import argparse
import pandas as pd
import bte
import numpy as np
import os
from github import Github
import glob

def argparser():
    parser = argparse.ArgumentParser(description="Write markdown issues representing the top N proposed sublineages of a report.")
    parser.add_argument("-i", "--input", required=True, help='Path to the lineage report (generate_lineage_report.py) to generate markdown issues from.')
    parser.add_argument("-t", "--tree", required=True, help='Path to the tree-containing protobuf the lineage report was generated from.')
    parser.add_argument("-m", "--metadata", required=True, help='Path to the metadata file corresponding to the tree the report was generate from.')
    parser.add_argument('-n', "--number", default=3, type=int, help="Write markdown reports for the top n lineages.")
    parser.add_argument("-s", "--sort", default='proposed_sublineage_score', help="Choose a column to sort by.")
    parser.add_argument('-p', "--prefix", default='proposal_', help="String to use as prefixes for output files.")
    parser.add_argument("-j", "--jsonsize", default=500, type=int, help="Maximum size of the json output for each sublineage.")
    parser.add_argument("-l", "--local", action='store_true', help="Set to write issues to local files only. Default posts issues to https://github.com/jmcbroome/auto-pango-designation/issues")
    parser.add_argument("-c", "--samplecount", default=15, type=int, help="Include up to this many sample names in the text of the report proposing the lineage.")
    parser.add_argument("-k", "--skip", action='store_true', help="Use to skip reporting any lineages that overlap with sample sets enumerated in *_samples.txt files in this directory.")
    args = parser.parse_args()
    return args
 
def get_date(d):
    try:
        return dt.datetime.strptime(d,"%Y-%m-%d")
    except:
        return np.nan

def write_report(row, prefix, samplenames, samplecount):
    fstr = []
    fstr.append("{} is a proposed child lineage of {} including {} samples.".format(row.proposed_sublineage, row.parent, row.proposed_sublineage_size))
    if type(row.earliest_child) != float and type(row.latest_child) != float:
        fstr.append("The earliest sample is {} and was found on {}, and the latest is {} and was found on {}.".format(samplenames[0], row.earliest_child, samplenames[1], row.latest_child))
    else:
        fstr.append("Dates could not be identified from the metadata for these samples.")
    if type(row.child_regions) != float and row.child_regions != np.nan: 
        ccount = row.child_regions.count(",") + 1
        common = row.child_regions.split(",")[0]
        commonprop = round(float(row.child_region_percents.split(",")[0])*100,2)
        if ccount == 1:
            fstr.append("It is found in {} only.".format(row.child_regions))
        else:
            fstr.append("It is found in {} countries, most commonly {} where {:.2f}% of its samples were sequenced.".format(ccount, common, commonprop))
    else:
        fstr.append("Countries could not be identified for these samples.")
    if row.net_escape_gain > 0:
        fstr.append("{} has a Bloom Lab escape score of {:.2f}, -{:.2f} over the parent lineage.".format(row.proposed_sublineage, row.sublineage_escape, row.net_escape_gain))
    else:
        fstr.append("{} has a Bloom Lab escape score of {:.2f}, unchanged from the parent lineage.".format(row.proposed_sublineage, row.sublineage_escape))
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
    fstr.append("The following samples are included: ")
    for s in samplenames:
        fstr.append(s)
    remainder = samplecount - len(samplenames)
    if remainder > 0:
        fstr.append("As well as {} additional samples.".format(remainder))
    return fstr

def write_sample_list(t, mdf, nid, name, prefix, count, skipset, use_skip):
    selection = []
    with open(prefix + name + "_samples.txt","w+") as outf:
        samples = t.get_leaves_ids(nid)
        if use_skip:
            if any([s for s in samples if s in skipset]):
                print("Proposal {} overlaps with existing proposals; skipping")
                return [], len(samples)
        metadata = mdf[mdf.strain.isin(samples)]
        if metadata.shape[0] == 0:
            print("None of the samples of lineage {} have accompanying metadata; skipping".format(name))
            return [], len(samples)
        # print(metadata.columns)
        metadata.apply(lambda row: print("\t".join([str(v) for v in [row.strain, row.country, row.date]],),file=outf),axis=1)
        selection.append(metadata.iloc[0].strain)
        selection.append(metadata.iloc[-1].strain)
        target = min(count, metadata.shape[0]-2)
        if target > 0:
            selection.extend(metadata.iloc[1:-1,:].strain.sample(target,replace=False))
    return selection, len(samples)

def get_current_proposed_covered():
    samset = set()
    sfiles = glob.glob("*_samples.txt")
    for f in sfiles:
        with open(f) as inf:
            for entry in inf:
                sample, country, date = entry.strip().split()
                samset.add(sample)
    return samset

def write_json(t, nid, parent_nid, name, prefix, size, metafile = None):
    outn = prefix + name + ".json"
    samples_to_use = t.get_leaves_ids(nid)
    parent_samples = t.get_leaves_ids(parent_nid)
    if len(samples_to_use) + len(parent_samples) > size:
        target = size - len(samples_to_use)
        parents_to_use = list(np.random.choice(parent_samples, size = target, replace = False))
    else:
        parents_to_use = parent_samples
    total = samples_to_use + parents_to_use
    if metafile != None:
        t.write_json(outn, samples = total, title = name, metafiles = [metafile])
    else:
        t.write_json(outn, samples = total, title = name, metafiles = [])

def main():
    args = argparser()
    t = bte.MATree(args.tree)
    tn = args.tree.split(".")[0]
    df = pd.read_csv(args.input,sep='\t')
    mdf = pd.read_csv(args.metadata,sep='\t')
    mdf['date'] = mdf.date.apply(get_date)
    mdf.sort_values('date',inplace=True)
    samset = get_current_proposed_covered()
    if len(samset) == 0:
        print("Found no current samples in *_samples.txt files.")
    else:
        print("Found {} samples in *_samples.txt files; fresh proposals containing any of these samples will not be reported as new lineages.".format(len(samset)))
    i = 0
    for ind,d in df.sort_values(args.sort).iterrows():
        if i >= args.number:
            break
        print("Recording output for {}".format(d.proposed_sublineage))
        print("Writing samples...")
        selected_samples, scount = write_sample_list(t, mdf, d.proposed_sublineage_nid, d.proposed_sublineage, args.prefix, args.samplecount, samset, args.skip)
        if len(selected_samples) == 0:
            print("{} has no dates or includes samples already covered by open proposals; skipping.".format(d.proposed_sublineage))
            continue
        print("Writing report...")
        report = write_report(d, args.prefix, selected_samples, scount)
        if args.local:
            with open(args.prefix + d.proposed_sublineage + ".md","w+") as outf:
                print("\n".join(report),file=outf)
        else:
            g = Github(os.getenv("API_KEY"))
            r = g.get_user().get_repo("auto-pango-designation")
            titlestring = "Sublineage {} of {}".format(d.proposed_sublineage, d.parent)
            r.create_issue(title=titlestring,body="\n".join(report))
        # print("Writing json...")
        # write_json(t, d.proposed_sublineage_nid, d.parent_nid, d.proposed_sublineage, args.prefix, args.jsonsize, args.metadata)
        i += 1
    with open(tn+".issues.log","w+") as outf:
        print("Produced files for {} lineage proposals.".format(args.number),file=outf)

if __name__ == "__main__":
    main()