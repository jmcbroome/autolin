import sys
sys.path.append("./SARS2_RBD_Ab_escape_maps/")
import bindingcalculator as bc
import bte
import pandas as pd
import numpy as np
import datetime as dt
import argparse
from urllib import parse

def argparser():
    parser = argparse.ArgumentParser(description="Compute detailed lineage reports for all existing lineages in the tree.")
    parser.add_argument("-i", "--input", required=True, help='Path to protobuf to compute reports from.')
    parser.add_argument("-p", "--proposed", required=True, help='Path to the file containing dumped sublineage proposals.')
    parser.add_argument("-o", "--output", help='Name of the output table.',default=None,required=True)
    parser.add_argument("-m", "--metadata", help="Path to a metadata file matching the protobuf.",required=True)
    parser.add_argument("-f", "--reference", default=None, help="Path to a reference fasta file. Use with -g to annotate amino acid changes and immune escape in the expanded output.")
    parser.add_argument("-g", "--gtf", default=None, help="Path to a reference gtf file. Use with -f to annotate amino acid changes and immune escape in the expanded output.")
    args = parser.parse_args()
    return args

def get_date(d):
    try:
        return dt.datetime.strptime(d,"%Y-%m-%d")
    except:
        return np.nan

def write_taxonium_url(parentlin, mutations):
    urlbase = 'https://taxonium.org/?backend=https://api.cov2tree.org&'
    searchbase = {"key":"aa1","type":"boolean","method":"boolean","text":parentlin,"gene":"S","position":484,"new_residue":"any","min_tips":0}
    searchbase['subspecs'] = [{"key":"ab0","type":"meta_pango_lineage_usher","method":"text_exact","text":parentlin,"gene":"S","position":484,"new_residue":"any","min_tips":0}]
    #taxonium uses key values that are distinct for every search but seemingly arbitrary.
    keyi = 1
    keys_used = []
    for gm in mutations:
        if ":" not in gm:
            print("WARNING: mutation {} failing to parse".format(gm))
            continue
        gene, m = gm.split(":")
        loc = m[1:-1]
        state = m[-1]
        key = str(keyi) + 'ab'
        keys_used.append(key)
        keyi += 1
        searchbase['subspecs'].append({"key":key,"type":"genotype","method":"genotype","text":"","gene":gene,"position":loc,"new_residue":state,"min_tips":0})
    searchbase['boolean_method'] = 'and'
    queries = parse.urlencode([("srch",'[' + "".join(str(searchbase).split()).replace("'",'"') + ']'),("enabled",'{"aa1":"true"}'),("zoomToSearch",0)])
    return urlbase + queries

def fill_output_table(t,pdf,mdf,fa_file=None,gtf_file=None):
    mdf.set_index('strain',inplace=True)
    def parent_lineage_size(row):
        samples = t.get_leaves_ids(row.parent_nid)
        return len(samples)
    def sublineage_size(row):
        subsamples = t.get_leaves_ids(row.proposed_sublineage_nid)
        return len(subsamples)
    print("Computing sublineage percentages.")
    pdf['parent_lineage_size'] = pdf.apply(parent_lineage_size,axis=1)
    pdf['proposed_sublineage_percent'] = round(pdf.proposed_sublineage_size/pdf.parent_lineage_size,2)
    def parsimony_parent(row):
        parent_parsimony = sum([len(n.mutations) for n in t.depth_first_expansion(row.parent_nid)])
        return parent_parsimony
    def parsimony_child(row):
        child_parsimony = sum([len(n.mutations) for n in t.depth_first_expansion(row.proposed_sublineage_nid)])
        return child_parsimony
    print("Computing parsimony percentages.")
    pdf['parent_parsimony'] = pdf.apply(parsimony_parent,axis=1)
    pdf['proposed_sublineage_parsimony'] = pdf.apply(parsimony_child,axis=1)
    pdf['parsimony_percent'] = round(pdf.proposed_sublineage_parsimony/pdf.parent_parsimony,2)
    mdf['date'] = mdf.date.apply(get_date)
    def get_start_ends(row):
        parent_samples = set(t.get_leaves_ids(row.parent_nid))
        child_samples = set(t.get_leaves_ids(row.proposed_sublineage_nid))
        parent_only = parent_samples - child_samples
        try:
            parent_dates = mdf.loc[[p for p in parent_only if p in mdf.index]].date
            child_dates = mdf.loc[[c for c in child_samples if c in mdf.index]].date
            return min(parent_dates),max(parent_dates),min(child_dates),max(child_dates)
        except KeyboardInterrupt:
            raise KeyboardInterrupt
        except IndexError:
            return np.nan,np.nan,np.nan,np.nan
    print("Computing start and end dates.")
    applied_pdf = pdf.apply(lambda row: get_start_ends(row), axis='columns', result_type='expand')
    pdf = pd.concat([pdf, applied_pdf], axis='columns')
    pdf = pdf.rename({0:'earliest_parent',1:'latest_parent',2:'earliest_child',3:'latest_child'},axis=1)
    pdf['log_score'] = np.log10(pdf.proposed_sublineage_score)
    print("Tracking country composition.")
    def get_regions(nid):
        try:
            return ",".join(list(mdf.loc[[l for l in t.get_leaves_ids(nid) if l in mdf.index]].country.value_counts().index))
        except:
            return np.nan
    def get_regions_percents(nid):        
        try:
            return ",".join([str(round(p,2)) for p in mdf.loc[[l for l in t.get_leaves_ids(nid) if l in mdf.index]].country.value_counts(normalize=True)])
        except:
            return np.nan
    pdf['child_regions'] = pdf.proposed_sublineage_nid.apply(get_regions)
    # pdf['ParentRegions'] = pdf.parent_nid.apply(get_regions)
    pdf['child_region_percents'] = pdf.proposed_sublineage_nid.apply(get_regions_percents)
    # pdf['ParentRegionPercents'] = pdf.parent_nid.apply(get_regions_percents)
    def host_jump(row):
        try:
            return mdf.loc[[l for l in t.get_leaves_ids(proposed_sublineage_nid) if l in mdf.index]].host.nunique() > 1
        except:
            return False
    print("Identifying host jumps.")
    pdf['host_jump'] = pdf.apply(host_jump,axis=1)
    print("Generating cov-spectrum URLs.")
    def generate_url(row):
        child_mset = t.get_haplotype(row.proposed_sublineage_nid)
        parent_mset = t.get_haplotype(row.parent_nid)
        net_mset = child_mset - parent_mset
        mset_str = "[{}-of:{}]".format(len(net_mset), ", ".join(list(net_mset)))
        query = parse.urlencode([('variantQuery','nextcladePangoLineage:' + row.parent + "*&" + mset_str)])
        url = "https://cov-spectrum.org/explore/World/AllSamples/AllTimes/variants?" + query
        return url
    pdf['link'] = pdf.apply(generate_url,axis=1)
    print("Collecting mutations.")
    def get_separating_mutations(row):
        hapstring = []
        for n in t.rsearch(row.proposed_sublineage_nid,True):
            if n.id == row.parent_nid:
                break
            hapstring.append(",".join(n.mutations))
        return ">".join(hapstring[::-1])
    pdf['mutations'] = pdf.apply(get_separating_mutations,axis=1)
    def get_growth_score(row):
        try:
            time = (row.latest_child - row.earliest_child).weeks + 1
            return np.sqrt(row.proposed_sublineage_size) / time
        except AttributeError:
            return np.nan
    pdf['growth_score'] = pdf.apply(get_growth_score,axis=1)
    if gtf_file != None and fa_file != None:
        print("Performing translation and computing antibody binding scores.")
        translation = t.translate(fasta_file = fa_file, gtf_file = gtf_file)
        calculator = bc.BindingCalculator(csv_or_url='SARS2_RBD_Ab_escape_maps/processed_data/escape_calculator_data.csv')
        hstrs = []
        cev = []
        pev = []
        nev = []
        for i, row in pdf.iterrows():
            hapstring = []
            child_changes = []
            parent_changes = []
            past_parent = False
            for n in t.rsearch(row.proposed_sublineage_nid,True):
                aas = [a for a in translation.get(n.id,[]) if a.original_aa != a.alternative_aa] #ignore synonymous changes.
                if n.id == row.parent_nid:
                    past_parent = True
                if past_parent:
                    parent_changes.extend([a.aa_index for a in aas if a.gene == 'S' and a.aa_index in calculator.sites])
                else:
                    child_changes.extend([a.aa_index for a in aas if a.gene == 'S' and a.aa_index in calculator.sites])
                    hapstring.append(",".join([aav.gene+":"+aav.aa for aav in aas]))
            hstr = ">".join(hapstring[::-1])
            child_escape = calculator.binding_retained(child_changes + parent_changes)
            parent_escape = calculator.binding_retained(parent_changes)
            net_escape_gain = parent_escape - child_escape
            hstrs.append(hstr)
            cev.append(child_escape)
            pev.append(parent_escape)
            nev.append(net_escape_gain)
        pdf['aa_changes'] = hstrs
        pdf['sublineage_escape'] = cev
        pdf['parent_escape'] = pev
        pdf['net_escape_gain'] = nev   
        def changes_to_list(aacstr):
            changes = []
            for n in aacstr.split(">"):
                if len(n) > 0:
                    changes.extend(n.split(","))
            return changes
        pdf['taxlink'] = pdf.apply(lambda row:write_taxonium_url(row.parent, changes_to_list(row.aa_changes)),axis=1)
    return pdf

def main():
    args = argparser()
    mdf = pd.read_csv(args.metadata,sep='\t')
    t = bte.MATree(args.input)
    pdf = pd.read_csv(args.proposed,sep='\t')
    odf = fill_output_table(t,pdf,mdf,args.reference,args.gtf)
    odf.to_csv(args.output,sep='\t',index=False)

if __name__ == "__main__":
    main()