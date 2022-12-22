import sys
sys.path.append("~/bin:")
import argparse
import pandas as pd
import bte
import numpy as np
import os
import datetime as dt
from github import Github
from pango_aliasor.aliasor import Aliasor
global_aliasor = Aliasor()

def argparser():
    parser = argparse.ArgumentParser(description="Open a pull request containing all lineages of a report added to a pango designation repository.")
    parser.add_argument("-r", "--repository", help='Path to the designation repository to be updated.', default = "./auto-pango-designation")
    parser.add_argument("-i", "--input", required=True, help='Path to the lineage report (from generate_lineage_report.py).')
    parser.add_argument("-t", "--tree", required=True, help='Path to the tree-containing protobuf the lineage report was generated from.')
    parser.add_argument("-l", "--local", action='store_true', help="Set to write issues to local files only. Default opens a pull request to https://github.com/jmcbroome/auto-pango-designation/issues from a custom branch")
    parser.add_argument("-c", "--representative", default=5000, type=int, help="Include up to this many representative samples for each lineage in lineages.csv.")
    parser.add_argument("-s", "--samples", default="None", help="Use only samples from the indicated file to update lineages.csv. Default behavior uses any samples.")
    parser.add_argument("-m", "--maximum",type=int,default=20,help="Include up to this many lineages in the body of the pull request. Default 20")
    parser.add_argument("-u", "--countries",type=int,default=1,help="Set to a minimum number of countries the lineage exists in to report it. Default is 1.")
    parser.add_argument("-o", "--output_report",default=None,help="Save the output report table to a tsv.")
    parser.add_argument("-a", "--active_since",default=None,help="Only report lineages actively sampled since the indicated date (formatted YYYY-MM-DD). Default reports all.")
    parser.add_argument("-d", "--samples_different",type=int,default=5,help="Ignore proposals which aren't at least this many samples smaller than the parent to prevent mostly overlapping lineage proposals.")
    parser.add_argument("--automerge", action='store_true', help='Immediately merge this pull request if permissions allow.')
    parser.add_argument("-e", "--metadata", default=None,help="Path to a sample-level metadata file. Required for -M")
    parser.add_argument("-M", "--model_growth", action='store_true', help='Infer an estimate for the exponential growth coefficient for all lineage proposals which pass other filters and contain sufficient geospatial diversity for inference.')
    parser.add_argument("--draws",type=int,default=1000,help="Set the draws parameter for the Bayesian growth model.")
    parser.add_argument("--tune", type=int,default=1500,help="Set the tuning parameter for the sampler.")
    parser.add_argument("--min_country_weeks",type=int,default=5,help="Require at least this many valid data points for inference.")
    parser.add_argument("--target_accept",type=float,default=.8,help="Set the target acceptance parameter for the sampler.")
    args = parser.parse_args()
    return args

def compress_lineage(al):
    #handling to deal with recompression of auto-extended names that may or may not be in the keys.
    fields = al.split(".")
    outer_fields = []
    while len(fields) > 0:
        name = ".".join(fields)
        try:
            fn = global_aliasor.compress(name) 
            if len(outer_fields) > 0:
                return fn + '.' + ".".join(reversed(outer_fields))
            else:
                return fn
        except KeyError:
            outer_fields.append(fields[-1])
            fields = fields[:-1]
    return al

def get_region_summary(row):
    regions_to_report = []
    cumprop = 0
    for cr, propstr in zip(row.child_regions.split(","),row.child_region_percents.split(",")):
        prop = float(propstr)
        regions_to_report.append(cr)
        cumprop += prop
        if cumprop > .8 or len(regions_to_report) >= 3:
            break
    if len(regions_to_report) == 1:
        cstr = regions_to_report[0]
    elif len(regions_to_report) == 2:
        cstr = regions_to_report[0] + " and " + regions_to_report[1]
    elif len(regions_to_report) == 3:
        cstr = ", ".join(regions_to_report[:-1]) + ", and " + regions_to_report[-1]
    else:
        cstr = ", ".join(regions_to_report[:3]) + f", and {len(regions_to_report)-3} additional countries."
    return cstr

def get_aa_set(row):
    aastr = []
    for aav in row.aa_path.split(">"):
        if len(aav) > 0:
            for aa in aav.split(","):
                al = aa.split(":")[1]
                oal = al[-1] + al[1:-1] + al[0]
                opp = aa.split(":")[0] + oal
                if opp not in aastr:
                    aastr.append(aa)
                else:
                    aastr.remove(opp)
    return aastr

def write_note(row):
    unalias = global_aliasor.uncompress(row.proposed_sublineage[5:])
    cstr = get_region_summary(row)
    aastr = get_aa_set(row)
    outstr = ['auto.' + compress_lineage(unalias) + "\t", "Alias of auto." + unalias]
    if len(aastr) > 0:
        outstr.append(", defined by " + ", ".join(aastr))
    if len(cstr) > 0:
        outstr.append(", found in " + cstr)
    outstr.append(". Automatically inferred by https://github.com/jmcbroome/automate-lineages-prototype.")
    return ''.join(outstr)

def get_reps(nid, t, target = 5000, allowed = set()):
    total = t.get_leaves_ids(nid)
    if len(allowed) > 0:
        total = [l for l in total if l in allowed]
    if len(total) <= target:
        return total
    else:
        return list(np.random.choice(total, replace=False, size=target))

def open_pr(branchname,trepo,automerge,reqname,pdf):
    g = Github(os.getenv("API_KEY"))
    repo = g.get_user().get_repo("auto-pango-designation")
    sb = repo.get_branch("master")
    if branchname not in [b.name for b in repo.get_branches()]:
        repo.create_git_ref(ref='refs/heads/' + branchname, sha=sb.commit.sha)
    for git_file in ["lineages.csv","lineage_notes.txt"]:
        print(f"Updating {git_file}...",file=sys.stderr)
        contents = repo.get_contents(git_file)
        with open(trepo+"/"+git_file) as inf:
            newcontent = inf.read()
        repo.update_file(contents.path, "Updating with new lineages.", newcontent, contents.sha, branch=branchname)
    repo.create_pull(title="New Lineages Update: " + reqname, body=pdf.to_markdown(index=False), head=branchname, base="master")
    if automerge:
        #merge to master, then delete this branch.
        head = repo.get_branch(branchname)
        repo.merge("master", head.commit.sha, "automerge")
        ref = repo.get_git_ref(f"heads/{branchname}")
        ref.delete()

def update_lineage_files(pdf, t, repo, rep, allowed, annotes):
    lincsv = repo + "/lineages.csv"
    skip = set()
    with open(lincsv, "a") as outf:
        for i,row in pdf.iterrows():
            sn = annotes.get(row.proposed_sublineage, None)
            if sn == None:
                print(f"WARNING: lineage {row.proposed_sublineage} not found on the input tree! Skipping")
                skip.add(row.proposed_sublineage)
                continue
            rsamples = get_reps(sn, t, rep, allowed)
            if len(rsamples) == 0:
                print(f"WARNING: no representative samples found for lineage {row.proposed_sublineage}! Skipping",file=sys.stderr)
                skip.add(row.proposed_sublineage)
                continue
            for rs in rsamples:
                print(rs + "," + row.proposed_sublineage, file=outf)
    print(f"{pdf.shape[0]-len(skip)} lineages added to lineages.csv; {len(skip)} skipped for having no high quality descendents.")
    notecsv = repo + "/lineage_notes.txt"
    pdf = pdf[~pdf.proposed_sublineage.isin(skip)]
    with open(notecsv, "a") as outf:
        pdf.apply(lambda row: print(write_note(row), file=outf), axis=1)
    print(f"Updated lineages.txt and lineages.csv with {pdf.shape[0]} additional lineages.")
    return pdf

def get_date(d):
    try:
        return dt.datetime.strptime(d,"%Y-%m-%d")
    except:
        return np.nan

def main():
    args = argparser()
    pdf = pd.read_csv(args.input,sep='\t')
    pdf = pdf[(pdf.child_regions_count >= args.countries) & (pdf.parent_lineage_size - pdf.proposed_sublineage_size >= args.samples_different)]
    if args.active_since != None and args.active_since != "None":
        pdf = pdf[(pdf.latest_child.apply(get_date) >= dt.datetime.strptime(args.active_since,"%Y-%m-%d"))]
    if args.model_growth:
        print("Modeling growth with PyMC3.",file=sys.stderr)
        from model_growth import get_growth_model
        if args.metadata == None:
            print("ERROR: -e must be set with -M output.")
            sys.exit(1)
        mdf = pd.read_csv(args.metadata,sep='\t')
        mdf = mdf[mdf.auto_annotation.isin(pdf.proposed_sublineage)]
        print(f"{mdf.shape[0]} samples to be used among {pdf.proposed_sublineage.nunique()} lineage models.")
        mdf['date'] = mdf.date.apply(lambda x:get_date(x))
        growd = get_growth_model(mdf, args.min_country_weeks, args.target_accept, args.tune, args.draws)
        pdf['Exponential Growth Coefficient CI'] = pdf.proposed_sublineage.apply(lambda x:growd.get(x,(np.nan,np.nan))) 
        pdf['Minimum Growth'] = pdf['Exponential Growth Coefficient CI'].apply(lambda x:x[0])
        pdf = pdf.sort_values("Minimum Growth",ascending=False)
    pdf = pdf.head(args.maximum)
    allowed = set()
    if args.samples != "None" and args.samples != None:
        with open(args.samples) as inf:
            for entry in inf:
                allowed.add(entry.strip())
        print(f"{len(allowed)} total samples are available for updating lineages.csv",file=sys.stderr)
    t = bte.MATree(args.tree)
    tannotes = t.dump_annotations()
    pdf = update_lineage_files(pdf, t, args.repository, args.representative, allowed, tannotes)
    pdf['link'] = pdf.link.apply(lambda x:f"[View On Cov-Spectrum]({x})")
    pdf['taxlink'] = pdf.taxlink.apply(lambda x:f"[View On Taxonium (Public Samples Only)]({x})")
    pdf['Regions'] = pdf.apply(get_region_summary,axis=1)
    def get_mutation_set(row):
        childhap = t.get_haplotype(row.proposed_sublineage_nid)
        parenthap = t.get_haplotype(row.parent_nid)
        return ",".join(list(childhap.difference(parenthap)))
    pdf['Nucleotide Changes'] = pdf.apply(get_mutation_set,axis=1)
    pdf['Amino Acid Changes'] = pdf.apply(lambda row: ",".join(get_aa_set(row)),axis=1)
    targets = ['proposed_sublineage', 'parent', 'proposed_sublineage_size','earliest_child','latest_child','Regions','Nucleotide Changes','Amino Acid Changes','link','taxlink']
    if "Exponential Growth Coefficient CI" in pdf.columns:
        targets.insert(3,'Exponential Growth Coefficient CI')
    pdf = pdf[targets]
    pdf = pdf.rename({"proposed_sublineage":"Lineage Name", "parent":"Parent Lineage", "proposed_sublineage_size":"Size","earliest_child":"Earliest Appearance","latest_child":"Latest Appearance","final_date":"Last Checked","child_regions":"Circulating In","link":"View On Cov-Spectrum","taxlink":"View On Taxonium (Public Samples Only)"},axis=1)
    if args.output_report != None:
        pdf.to_csv(args.output_report,index=False,sep='\t')
    if not args.local:
        branch = str(pdf.shape[0]) + "_" + "_".join(args.tree.split("/")[-1].split(".")[:2])
        date = args.tree.split("/")[-1].split(".")[1]
        reqname = f"{str(pdf.shape[0])} new lineages active through {date}"
        open_pr(branch, args.repository, args.automerge, reqname, pdf)
        print("Github updated.")
if __name__ == "__main__":
    main()