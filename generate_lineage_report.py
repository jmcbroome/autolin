import bte
import pandas as pd
import numpy as np
import datetime as dt
import argparse

def argparser():
    parser = argparse.ArgumentParser(description="Compute detailed lineage reports for all existing lineages in the tree.")
    parser.add_argument("-i", "--input", required=True, help='Path to protobuf to compute reports from.')
    parser.add_argument("-p", "--proposed", required=True, help='Path to the file containing dumped sublineage proposals.')
    parser.add_argument("-o", "--output", help='Name of the output table.',default=None,required=True)
    parser.add_argument("-m", "--metadata", help="Path to a metadata file matching the protobuf.",required=True)
    args = parser.parse_args()
    return args

def get_date(d):
    try:
        return dt.datetime.strptime(d,"%Y-%m-%d")
    except:
        return np.nan

def is_successive(row):
    if row.EarliestChild > row.EarliestParent and row.LatestChild > row.LatestParent:
        return True
    else:
        return False

def fill_output_table(t,pdf,mdf):
    mdf.set_index('strain',inplace=True)
    def sublineage_size(row):
        samples = t.get_leaves_ids(row.parent_nid)
        subsamples = t.get_leaves_ids(row.proposed_sublineage_nid)
        return len(subsamples)/len(samples)
    print("Computing sublineage percent")
    pdf['SublineagePercent'] = pdf.apply(sublineage_size,axis=1)
    def parsimony_percent(row):
        parent_parsimony = sum([len(n.mutations) for n in t.depth_first_expansion(row.parent_nid)])
        child_parsimony = sum([len(n.mutations) for n in t.depth_first_expansion(row.proposed_sublineage_nid)])
        if parent_parsimony == 0:
            return np.nan
        else:
            return child_parsimony/parent_parsimony
    print("Computing parsimony percent")
    pdf['ParsimonyPercent'] = pdf.apply(parsimony_percent,axis=1)
    def get_start_ends(row):
        parent_samples = set(t.get_leaves_ids(row.parent_nid))
        child_samples = set(t.get_leaves_ids(row.proposed_sublineage_nid))
        parent_only = parent_samples - child_samples
        try:
            parent_dates = mdf.loc[list(parent_only)].date
            child_dates = mdf.loc[list(child_samples)].date
            return min(parent_dates),max(parent_dates),min(child_dates),max(child_dates)
        except KeyboardInterrupt:
            raise KeyboardInterrupt
        except:
            return np.nan,np.nan,np.nan,np.nan
    print("Computing start and end dates")
    applied_pdf = pdf.apply(lambda row: get_start_ends(row), axis='columns', result_type='expand')
    pdf = pd.concat([pdf, applied_pdf], axis='columns')
    pdf = pdf.rename({0:'EarliestParent',1:'LatestParent',2:'EarliestChild',3:'LatestChild'},axis=1)
    pdf['LogScore'] = np.log10(pdf.proposed_sublineage_score)
    pdf["Successive"] = pdf.apply(is_successive,axis=1)
    mdf['date'] = mdf.date.apply(get_date)
    print("Checking geography")
    def is_founder(row):
        try:
            parent_regions = mdf.loc[t.get_leaves_ids(row.parent_nid)].country.value_counts()
            child_regions = mdf.loc[t.get_leaves_ids(row.proposed_sublineage_nid)].country.value_counts()
            return child_regions.index[0] in parent_regions.index
        except:
            return False
    pdf['IsFounder'] = pdf.apply(is_founder,axis=1)
    print("Doing international")
    def is_international(row):
        try:
            return mdf.loc[t.get_leaves_ids(row.proposed_sublineage_nid)].country.nunique() > 1
        except:
            return False
    pdf['International'] = pdf.apply(is_international, axis=1)
    def generate_url(nid):
        mset = t.get_haplotype(nid)
        url = "https://cov-spectrum.org/explore/World/AllSamples/AllTimes/variants?variantQuery=[" + str(len(mset)) + "-of:"
        start = True
        for m in mset:
            if start:
                start = False
            else:
                url += ', '
            url += m
        url += ']'
        return url
    pdf['Links'] = pdf.proposed_sublineage_nid.apply(generate_url)
    return pdf

def main():
    args = argparser()
    mdf = pd.read_csv(args.metadata,sep='\t')
    t = bte.MATree(args.input)
    pdf = pd.read_csv(args.proposed,sep='\t')
    odf = fill_output_table(t,pdf,mdf)
    odf.to_csv(args.output,sep='\t',index=False)

if __name__ == "__main__":
    main()