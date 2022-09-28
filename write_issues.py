import argparse
import pandas as pd
import bte

def argparser():
    parser = argparse.ArgumentParser(description="Write markdown issues representing the top N proposed sublineages of a report.")
    parser.add_argument("-i", "--input", required=True, help='Path to the lineage report (generate_lineage_report.py) to generate markdown issues from.')
    parser.add_argument("-t", "--tree", required=True, help='Path to the tree-containing protobuf the lineage report was generated from.')
    parser.add_argument("-m", "--metadata", required=True, help='Path to the metadata file corresponding to the tree the report was generate from.')
    parser.add_argument('-n', "--number", default=3, help="Write markdown reports for the top n lineages.")
    parser.add_argument("-s", "--sort", default='proposed_sublineage_score', help="Choose a column to sort by.")
    parser.add_argument('-p', "--prefix", default='proposal_', help="String to use as prefixes for output files.")
    parser.add_argument("-j", "--jsonsize", default=4000, help="Maximum size of the json output for each sublineage.")
    args = parser.parse_args()
    return args

def write_report(row, prefix):
    outf = open(prefix + row.proposed_sublineage + ".md","w+")
    print("{} is a proposed child lineage of {} including {} samples.".format(row.proposed_sublineage, row.parent, row.proposed_sublineage_size),file=outf)
    if row.earliest_child != np.nan and row.latest_child != np.nan:
        print("The earliest sample was found on {} and the latest on {}.".format(row.earliest_child, row.latest_child),file=outf)
    else:
        print("Dates could not be identified for these samples.")
    if row.child_regions != np.nan:
        ccount = row.child_regions.count(",") + 1
        common = row.child_regions.split(",")[0]
        commonprop = round(float(row.child_region_percents.split(",")[0])*100,2)
        if ccount == 1:
            print("It is found in {} only.".format(row.child_regions),file=outf)
        else:
            print("It is found in {} countries, most commonly {} where {} of its samples were sequenced.".format(ccount, common, commonprop),file=outf)
    else:
        print("Countries could not be identified for these samples.",file=outf)
    print("{} has a Bloom Lab escape score of {}, -{} over the parent lineage.".format(row.proposed_sublineage, row.sublineage_escape, row.net_escape_gain),file=outf)
    spikes = []
    total = 0
    for n in row.aa_changes.split(">"):
        for m in n.split(","):
            if m[0] == 'S':
                spikes.append(m)
            total += 1
    if len(spikes) == 0:
        print("It is not defined by any spike protein changes. It has {} defining protein changes overall.".format(total),file=outf)
    else:
        print("It is defined by the following spike protein changes: {}. There are {} defining protein changes overall.".format(",".join(spikes), total),file=outf)
    if row.host_jump:
        print("It represents a zoonotic event!",file=outf)
    print("View it on cov-spectrum: {}",file=outf)
    outf.close()

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
            print('\t'.join([str(v) for v in [s, country, date]]) file = outf)

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
        t.write_json(outn.encode("UTF-8"), samples = samples_to_use + parents_to_use, title = name.encode("UTF-8"), metafiles = [metafile])
    else:
        t.write_json(outn.encode("UTF-8"), samples = samples_to_use + parents_to_use, title = name.encode("UTF-8"), metafiles = [])

def main():
    args = argparser()
    t = bte.MATree(args.tree)
    tn = args.tree.split(".")[0]
    df = pd.read_csv(args.input,sep='\t')
    mdf = pd.read_csv(args.metadata,sep='\t').set_index('strain')
    for i,d in df.sort_values(args.sort)[:args.number]:
        write_report(d, args.prefix)
        write_sample_list(t, mdf, d.proposed_sublineage_nid, d.proposed_sublineage, args.prefix)
        submeta_name = d.proposed_sublineage + "_metadata.tsv"
        sdf = mdf.loc[t.get_leaves_ids(d.proposed_sublineage_nid)]
        if sdf.shape[0] > 0:
            sdf.reset_index().to_csv(submeta_name,sep='\t',index=False)
            write_json(t, d.proposed_sublineage_nid, d.parent_nid, d.proposed_sublineage, args.prefix, args.jsonsize, submeta_name)
        else:
            write_json(t, d.proposed_sublineage_nid, d.parent_nid, d.proposed_sublineage, args.prefix, args.jsonsize)
    with open(tn+"_issues.log","w+") as outf:
        print("Produced files for {} lineage proposals.".format(args.number),file=outf)

if __name__ == "__main__":
    main()