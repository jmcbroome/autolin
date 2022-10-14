import argparse
import pandas as pd
import bte
import numpy as np
import os
import sys
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
    parser.add_argument("--automerge", action='store_true', help='Immediately merge this pull request if permissions allow.')
    args = parser.parse_args()
    return args

def write_note(row):
    unalias = global_aliasor.uncompress(row.proposed_sublineage)
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
    else:
        cstr = ", ".join(regions_to_report[:-1]) + ", and " + regions_to_report[-1]
    aastr = []
    for aav in row.aa_changes.split(">"):
        if len(aav) > 0:
            for aa in aav.split(","):
                al = aa.split(":")[1]
                oal = al[-1] + al[1:-1] + al[0]
                opp = aa.split(":")[0] + oal
                if opp not in aastr:
                    aastr.append(aa)
                else:
                    aastr.remove(opp)
    outstr = ["Alias of " + unalias]
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
    # samples = []
    # for node in t.breadth_first_expansion(nid):
    #     if node.is_leaf():
    #         if len(allowed) == 0 or node.id in allowed:
    #             samples.append(node.id)
    #     if len(samples) >= target:
    #         break
    return samples

def open_pr(branchname,trepo,automerge):
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
    repo.create_pull(title="New Lineages Update: " + branchname, body="", head=branchname, base="master")
    if automerge:
        #merge to master, then delete this branch.
        head = repo.get_branch(branchname)
        repo.merge("master", head.commit.sha, "automerge")
        ref = repo.get_git_ref(f"heads/{branchname}")
        ref.delete()

def main():
    args = argparser()
    pdf = pd.read_csv(args.input,sep='\t')
    allowed = set()
    if args.samples != "None" and args.samples != None:
        with open(args.samples) as inf:
            for entry in inf:
                allowed.add(entry.strip())
        print(f"{len(allowed)} total samples are available for updating lineages.csv",file=sys.stderr)
    t = bte.MATree(args.tree)
    lincsv = args.repository + "/lineages.csv"
    skip = set()
    with open(lincsv, "a") as outf:
        for i,row in pdf.iterrows():
            sn = row.proposed_sublineage_nid
            rsamples = get_reps(sn, t, args.representative, allowed)
            if len(rsamples) == 0:
                print("WARNING: no representative samples found for lineage {}! Skipping".format(row.proposed_sublineage),file=sys.stderr)
                skip.add(row.proposed_sublineage)
                continue
            for rs in rsamples:
                print(rs + "," + row.proposed_sublineage, file=outf)
    print(f"{pdf.shape[0]-len(skip)} lineages added to lineages.csv; {len(skip)} skipped for having no high quality descendents.")
    notecsv = args.repository + "/lineage_notes.txt"
    pdf = pdf[~pdf.proposed_sublineage.isin(skip)]
    with open(notecsv, "a") as outf:
        pdf.apply(lambda row: print(write_note(row), file=outf), axis=1)
    print(f"Updated lineages.txt and lineages.csv with {pdf.shape[0]} additional lineages.")
    if not args.local:
        open_pr(str(pdf.shape[0]) + "_" + args.tree, args.repository, args.automerge)
        print("Github updated.")
if __name__ == "__main__":
    main()