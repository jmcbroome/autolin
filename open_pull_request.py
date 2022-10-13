import argparse
import pandas as pd
import bte
import numpy as np
import os
from github import Github
import sys
from pango_aliasor.aliasor import Aliasor
global_aliasor = Aliasor()

def argparser():
    parser = argparse.ArgumentParser(description="Open a pull request containing all lineages of a report added to a pango designation repository.")
    parser.add_argument("-r", "--repository", help='Path to the designation repository to be updated.', default = "./auto-pango-designation")
    parser.add_argument("-i", "--input", required=True, help='Path to the lineage report (from generate_lineage_report.py).')
    parser.add_argument("-t", "--tree", required=True, help='Path to the tree-containing protobuf the lineage report was generated from.')
    parser.add_argument('-p', "--prefix", default='proposal_', help="String to use as prefixes for output files for local output.")
    parser.add_argument("-l", "--local", action='store_true', help="Set to write issues to local files only. Default opens a pull request to https://github.com/jmcbroome/auto-pango-designation/issues from a custom branch")
    parser.add_argument("-c", "--representative", default=5000, type=int, help="Include up to this many representative samples for each lineage in lineages.csv.")
    parser.add_argument("-s", "--samples", default="None", help="Use only samples from the indicated file to update lineages.csv. Default behavior uses any samples.")
    parser.add_argument("--automerge", action='store_true', help='Immediately merge this pull request if permissions allow.')
    args = parser.parse_args()
    return args

def write_note(row):
    print(row,file=sys.stderr)
    unalias = global_aliasor.uncompress(row.proposed_sublineage)
    regions_to_report = []
    cumprop = 0
    for cr, propstr in zip(row.child_regions.split(","),row.child_region_percents.split(",")):
        prop = float(propstr)
        regions_to_report.append(cr)
        cumprop += prop
        if cumprop > .8 or len(regions_to_report) >= 5:
            break
    cstr = ",".join(regions_to_report[:-1]) + ", and " + regions_to_report[-1]
    aastr = []
    print(row.aa_changes,file=sys.stderr)
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
        outstr.append(", defined by " + ",".join(aastr))
    if len(cstr) > 0:
        outstr.append(", mostly from " + cstr)
    outstr.append(". Automatically inferred by https://github.com/jmcbroome/automate-lineages-prototype.")
    return ''.join(outstr)

def get_reps(nid, t, target = 5000, allowed = set()):
    samples = []
    for node in t.breadth_first_expansion(nid):
        if node.is_leaf():
            if len(allowed) == 0 or node.id in allowed:
                samples.append(node.id)
        if len(samples) >= target:
            break
    return samples

def open_pr(branchname,trepo,automerge):
    g = Github(os.getenv("API_KEY"))
    repo = g.get_user().get_repo("auto-pango-designation")
    sb = repo.get_branch("master")
    repo.create_git_ref(ref='refs/heads/' + branchname, sha=sb.commit.sha)
    for git_file in [trepo+"/lineages.csv",trepo+"/lineage_notes.txt"]:
        contents = repo.get_contents(git_file)
        with open(git_file) as inf:
            newcontent = inf.read()
        repo.update_file(contents.path, "Updating with new lineages.", newcontent, contents.sha, branch=branchname)
    repo.create_pull(title="New Lineages Update: " + branchname, body="", head=branchname, base="master")
    if automerge:
        head = repo.get_branch(branchname)
        repo.merge("master", head.commit.sha, "automerge")        

def main():
    args = argparser()
    pdf = pd.read_csv(args.input,sep='\t')
    print(pdf.columns,file=sys.stderr)
    notecsv = args.repository + "/lineage_notes.txt"
    with open(notecsv, "a") as outf:
        pdf.apply(lambda row: print(write_note(row), file=outf), axis=1)
    allowed = set()
    if args.samples != "None" and args.samples != None:
        with open(args.samples) as inf:
            for entry in inf:
                allowed.add(entry.strip())
        print(f"{len(allowed)} total samples are available for updating lineages.csv")
    t = bte.MATree(args.tree)
    lincsv = args.repository + "/lineages.csv"
    with open(lincsv, "a") as outf:
        for i,row in pdf.iterrows():
            sn = row.proposed_sublineage_nid
            rsamples = get_reps(sn, t, args.representative, allowed)
            for rs in rsamples:
                print(rs + "," + row.proposed_sublineage, file=outf)
    print(f"Updated lineages.txt and lineages.csv with {pdf.shape[0]} additional lineages.")
    if not args.local:
        open_pr(args.prefix + str(pdf.shape[0]), args.repository, args.automerge)
        print("Github updated.")
if __name__ == "__main__":
    main()