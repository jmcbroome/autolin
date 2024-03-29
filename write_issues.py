import sys
sys.path.append("~/bin:")
import argparse
import pandas as pd
import bte
import numpy as np
import os
from github import Github
import glob
import datetime as dt

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

def write_report(row, prefix, samplenames, samplecount, treename = None, treenode = None):
    fstr = []
    fstr.append("{} is a proposed sublineage of {} that includes {} samples.".format(row.proposed_sublineage, row.parent, row.proposed_sublineage_size))
    if treename != None and treenode != None:
        fstr.append("\nIt was identified from the tree file {} as the descendents of {}.".format(treename, treenode))
    if type(row.earliest_child) != float and type(row.latest_child) != float:
        fstr.append("\nThe earliest dated sample was found on {}.".format(row.earliest_child))
        fstr.append("\nThe latest dated sample was found on {}.".format(row.latest_child))
    else:
        fstr.append("\nDates could not be identified from the metadata for these samples.")
    common = None
    if type(row.child_regions) != float and row.child_regions != np.nan: 
        ccount = row.child_regions.count(",") + 1
        common = row.child_regions.split(",")[0]
        commonprop = round(float(row.child_region_percents.split(",")[0])*100,2)
        if ccount == 1:
            fstr.append("\nIt is found in {} only.".format(row.child_regions))
        else:
            fstr.append("\nIt is found in {} countries, most commonly {} where {:.2f}% of its samples were sequenced.".format(ccount, common, commonprop))
    else:
        fstr.append("\nCountries could not be identified for these samples.")
    if row.net_escape_gain > 0:
        fstr.append("\n{} has a Bloom Lab escape score of {:.2f}, -{:.2f} over the parent lineage.".format(row.proposed_sublineage, row.sublineage_escape, row.net_escape_gain))
    else:
        fstr.append("\n{} has a Bloom Lab escape score of {:.2f}, unchanged from the parent lineage.".format(row.proposed_sublineage, row.sublineage_escape))
    if type(row.mutations) != float:
        fstr.append("\nFull mutational path from the ancestral lineage to the proposed sublineage: {}".format(row.mutations))
    else:
        fstr.append("\nWARNING: No mutations are associated with this proposal. Check parameters!")
    if type(row.aa_changes) != float:
        fstr.append("\nFull amino acid change path from the ancestral lineage to the proposed sublineage: {}".format(row.aa_changes))
    else:
        fstr.append("\nNo amino acid changes occurred between the ancestral and the proposed sublineage.")
    # total_mutations = 0
    # last_mutations = []
    # if type(row.mutations) != float:
    #     iterate = row.mutations.split(">")
    #     lastchanges = iterate[-1].split(",")
    #     for n in iterate:
    #         for m in n.split(","):
    #             if len(m) > 0:
    #                 total_mutations += 1
    #     last_mutations = iterate[-1].split(",")
    # if len(last_mutations) == 0:
    #     fstr.append("\nNo mutations are associated directly with the lineage root node. It has {} mutations from the parent lineage.".format(total_mutations))
    # else:
    #     fstr.append("\nThe following mutations are directly associated with the lineage root node: {}. There are {} mutations from the parent lineage overall.".format(",".join(last_mutations), total_mutations))
    fortitle = []
    nonspike_fortitle_count = 0
    # total = 0
    # lastchanges = []
    if type(row.aa_changes) != float:
        iterate = row.aa_changes.split(">")
        # lastchanges = [i for i in iterate[-1].split(",") if i != ""]
        for n in iterate:
            for m in n.split(","):
                if len(m) > 0:
                    if m[0] == 'S':
                        fortitle.append(m)
                    elif nonspike_fortitle_count <= 3:
                        fortitle.append(m)
                        nonspike_fortitle_count += 1
                    # total += 1
    # if len(lastchanges) == 0:
    #     fstr.append("\nNo protein changes are associated directly with the lineage root node. It has {} protein changes from the parent lineage.".format(total))
    # else:
    #     fstr.append("\nThe following protein changes are directly associated with the lineage root node: {}. There are {} defining protein changes from the parent lineage overall.".format(",".join(lastchanges), total))
    if row.host_jump:
        fstr.append("\nIt represents a zoonotic event!")
    fstr.append("\nNOTE: The following links search by genotype and parent lineage; they may additionally highlight samples outside of the proposed clade that convergently evolved the same mutations or that were added after this lineage was inferred, ")
    fstr.append("or not display or highlight samples which reverted some of the key changes or were removed after this lineage was inferred. They may also take a moment to complete the search on loading.")
    fstr.append("\nView it on [cov-spectrum]({})".format(row.link))
    fstr.append("\nView it on [taxonium]({})".format(row.taxlink))
    fstr.append("\nThe following samples are included in this proposal: ")
    for s in samplenames:
        fstr.append(s)
    remainder = samplecount - len(samplenames)
    if remainder > 0:
        fstr.append("As well as {} additional samples.".format(remainder))
    return fstr, ','.join(fortitle), common

def write_sample_list(t, mdf, nid, name, prefix, count, skipset):
    selection = []
    with open(prefix + name + "_samples.txt","w+") as outf:
        # samples = t.get_leaves_ids(nid)
        #ensure that among the sample names printed, at least two have their MRCA at the node root. 
        samples = []
        childset = {}
        for child in t.get_node(nid).children:
            cs = t.get_leaves_ids(child.id)
            for s in cs:
                samples.append(s)
                childset[s] = child.id
            # samples.extend(cs)
            # childset[child.id] = set(cs)
        if any([s for s in samples if s in skipset]):
            print("Proposal {} overlaps with existing proposals; skipping")
            return [], len(samples)
        metadata = mdf[mdf.strain.isin(samples)]
        metadata['FromChild'] = metadata.strain.apply(lambda x:childset[s])
        if metadata.shape[0] == 0:
            print("None of the samples of lineage {} have accompanying metadata; skipping".format(name))
            return [], len(samples)
        def print_sample_info(row):
            try:
                datestr = row.date.strftime("%Y-%m-%d")
            except ValueError:
                datestr = "NaN"
            print("\t".join([str(v) for v in [row.strain, row.country, datestr]],),file=outf)
        metadata.apply(print_sample_info,axis=1)
        #pick one from each child group. This ensures that the LCA of the named group is the correct root of the lineage.
        selection.extend(metadata.groupby("FromChild", group_keys=False).strain.sample(1)) 
        remaining = metadata[~metadata.strain.isin(selection)]
        selection.extend(remaining.strain.sample(min([count - len(selection), remaining.shape[0]]), replace=False))
        assert len(selection) <= count
    return selection, len(samples)

def get_current_proposed_covered():
    samset = set()
    sfiles = glob.glob("*_samples.txt")
    for f in sfiles:
        with open(f) as inf:
            for entry in inf:
                try:
                    sample, country, date = entry.strip().split()
                    samset.add(sample)
                except:
                    print("WARNING: following did not parse correctly")
                    print(entry)
                    continue
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
    tn = ".".join(args.tree.split(".")[:-2])
    df = pd.read_csv(args.input,sep='\t')
    mdf = pd.read_csv(args.metadata,sep='\t')
    mdf['date'] = mdf.date.apply(get_date)
    mdf.sort_values('date',inplace=True)
    if args.skip:
        samset = get_current_proposed_covered()
    else:
        samset = set()
    if len(samset) == 0 and args.skip:
        print("Found no current samples in *_samples.txt files.")
    elif args.skip:
        print("Found {} samples in *_samples.txt files; fresh proposals containing any of these samples will not be reported as new lineages.".format(len(samset)))
    i = 0
    for ind,d in df.sort_values(args.sort, ascending=False).iterrows():
        if i >= args.number:
            break
        if "(" in d.proposed_sublineage:
            print("Skipping {} for appearing not to be from a pangolin lineage...".format(d.proposed_sublineage))
            continue
        print("Recording output for {}".format(d.proposed_sublineage))
        print("Writing samples...")
        selected_samples, scount = write_sample_list(t, mdf, d.proposed_sublineage_nid, d.proposed_sublineage, args.prefix, args.samplecount, samset)
        if len(selected_samples) == 0:
            print("{} has no metadata or includes samples already covered by open proposals; skipping.".format(d.proposed_sublineage))
            continue
        print("Writing report...")
        report, fortitle, common_country = write_report(d, args.prefix, selected_samples, scount, args.tree, d.proposed_sublineage_nid)
        titlestring = "{} {} samples".format(d.proposed_sublineage_size, d.parent)
        if common_country != None:
            titlestring += " in {}".format(common_country)
        if len(fortitle) > 0:
            titlestring += " with {}".format(fortitle)
        if args.local:
            with open(args.prefix + d.proposed_sublineage + ".md","w+") as outf:
                print("Title: " + titlestring, file=outf)
                print("\n".join(report), file=outf)
        else:
            g = Github(os.getenv("API_KEY"))
            r = g.get_user().get_repo("auto-pango-designation")
            r.create_issue(title=titlestring,body="\n".join(report))
        # print("Writing json...")
        # write_json(t, d.proposed_sublineage_nid, d.parent_nid, d.proposed_sublineage, args.prefix, args.jsonsize, args.metadata)
        i += 1
    with open(tn+".issues.log","w+") as outf:
        print("Produced files for {} lineage proposals.".format(args.number),file=outf)

if __name__ == "__main__":
    main()