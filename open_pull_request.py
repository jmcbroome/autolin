import argparse
import pandas as pd
import bte
import numpy as np
import os
from github import Github
from pango_aliasor.aliasor import Aliasor
global_aliasor = Aliasor()

def argparser():
    parser = argparse.ArgumentParser(description="Open a pull request containing all lineages of a report added to a pango designation repository.")
    parser.add_argument("-r", "--repository", help='Path to the designation repository to be updated.', default = "./auto-pango-designation")
    parser.add_argument("-i", "--input", required=True, help='Path to the lineage report (from generate_lineage_report.py).')
    parser.add_argument("-t", "--tree", required=True, help='Path to the tree-containing protobuf the lineage report was generated from.')
    parser.add_argument("-m", "--metadata", required=True, help='Path to the metadata file corresponding to the tree the report was generate from.')
    parser.add_argument('-p', "--prefix", default='proposal_', help="String to use as prefixes for output files for local output.")
    parser.add_argument("-l", "--local", action='store_true', help="Set to write issues to local files only. Default opens a pull request to https://github.com/jmcbroome/auto-pango-designation/issues from a custom branch")
    parser.add_argument("-c", "--representative", default=5000, type=int, help="Include up to this many representative samples for each lineage in lineages.csv.")
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
        if cumprop > .8 or len(regions_to_report) >= 5:
            break
    cstr = ",".join(regions_to_report[:-1]) + ", and " + regions_to_report[-1]
    aastr = []
    for aav in row.aa_changes.split(">"):
        for aa in aav.split(","):
            al = aa.split(":")[1]
            oal = al[-1] + al[1:-1] + al[0]
            opp = aa.split(":")[0] + oal
            if opp not in aastr:
                aastr.append(aa)
            else:
                aastr.remove(opp)
    outstr = "Alias of " + unalias
    if len(aastr) > 0:
        outstr.append(", defined by " + aastr)
    if len(cstr) > 0:
        outstr.append(", mostly from " + cstr)
    outstr.append(". Automatically inferred by https://github.com/jmcbroome/automate-lineages-prototype.")
    return outstr

def get_reps(nid, t, target = 5000):
    samples = []
    for node in t.breadth_first_expansion(nid):
        if node.is_leaf():
            samples.append(node.id)
        if len(samples) >= target:
            break
    return samples

def open_pr(branchname,automerge):
    g = Github(os.getenv("API_KEY"))
    repo = g.get_user().get_repo("auto-pango-designation")
    sb = repo.get_branch("master")
    repo.create_git_ref(ref='refs/heads/' + branchname, sha=sb.commit.sha)
    for git_file in ["lineages.csv","lineage_notes.txt"]:
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
    t = bte.MATree(args.tree)
    lincsv = args.repository + "/lineages.csv"
    notecsv = args.repository + "/lineage_notes.txt"
    pdf = pd.read_csv(args.input)
    with open(notecsv, "a") as outf:
        pdf.apply(lambda row: print(write_note(row), file=outf))
    with open(lincsv, "a") as outf:
        for sn in pdf.proposed_sublineage_nid:
            rsamples = get_reps(sn, t, args.representative)
            for rs in rsamples:
                print(rs + "," + row.proposed_sublineage, file=outf)
    # if not args.local:
        # print("ERROR: remote posting NYI")
        # exit(1)
    # else:
    #     with open(args.prefix + "prbody.txt","w+") as outf:
    #         print(write_pr_body(pdf), file=outf)
    open_pr(args.prefix + str(pdf.shape[0]), args.repository, args.automerge)
    print(f"Updates lineages.txt and lineages.csv with {pdf.shape[0]} additional lineages.")

if __name__ == "__main__":
    main()