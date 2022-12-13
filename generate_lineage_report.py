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
    for m in mutations:
        # if ":" not in gm:
            # print("WARNING: mutation {} failing to parse".format(gm))
            # continue
        # gene, m = gm.split(":")
        gene = 'nt' #use the nucleotide search to generate links, avoid inconsistency with different GTFs.
        loc = m[1:-1]
        state = m[-1]
        key = str(keyi) + 'ab'
        keys_used.append(key)
        keyi += 1
        searchbase['subspecs'].append({"key":key,"type":"genotype","method":"genotype","text":"","gene":gene,"position":loc,"new_residue":state,"min_tips":0})
    searchbase['boolean_method'] = 'and'
    queries = parse.urlencode([("srch",'[' + "".join(str(searchbase).split()).replace("'",'"') + ']'),("enabled",'{"aa1":"true"}'),("zoomToSearch",0)])
    return urlbase + queries

def compute_stratified_growth(mdf):
    mdf['date'] = mdf.date.apply(get_date)
    rc = mdf.groupby(['country','autolin',pd.Grouper(key='fdate', freq='1W')]).strain.count().reset_index()
    rc = rc.rename({"strain":"count"},axis=1).sort_values("fdate")
    rc['cumcount'] = rc.groupby(['country','autolin'])['count'].cumsum()
    rc['abs_strat_growth'] = rc.groupby(['country','autolin'])['cumcount'].diff()
    rc['strat_growth'] = rc.groupby(['country','autolin'])['cumcount'].pct_change()
    growdf = rc[(rc.abs_strat_growth > 5)].replace(np.inf, np.nan).dropna().groupby(['pango_lineage_usher','country']).strat_growth.describe()
    gp = growdf.groupby("autolin")
    return gp['mean'].mean().to_dict()

def fill_output_table(t,pdf,mdf,fa_file=None,gtf_file=None):
    print("Filling out metadata with terminal lineages.")
    def get_latest_lineage(s):
        for anc in t.rsearch(s):
            if len(anc.annotations) > 0:
                if len(anc.annotations[0]) > 0:
                    return anc.annotations[0]
    mdf['autolin'] = mdf.strain.apply(get_latest_lineage)
    print("Computing geographically stratified growth values.")
    lingrow = compute_stratified_growth(mdf)
    pdf['mean_stratified_growth'] = pdf.proposed_sublineage.apply(lambda x:lingrow.get(x,np.nan))
    mdf.set_index('strain',inplace=True)
    #parent lineage size has to be inclusive to get a sensible percentage.
    def parent_lineage_size(row):
        samples = t.get_leaves_ids(row.parent_nid)
        return len(samples)
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
    # def get_start_ends(row):
    #     parent_samples = set(t.get_leaves_ids(row.parent_nid))
    #     child_samples = set(t.get_leaves_ids(row.proposed_sublineage_nid))
    #     parent_only = parent_samples - child_samples
    #     try:
    #         parent_dates = mdf.loc[[p for p in parent_only if p in mdf.index]].date
    #         child_dates = mdf.loc[[c for c in child_samples if c in mdf.index]].date
    #         return min(parent_dates),max(parent_dates),min(child_dates),max(child_dates)
    #     except KeyboardInterrupt:
    #         raise KeyboardInterrupt
    #     except IndexError:
    #         return np.nan,np.nan,np.nan,np.nan
    def get_start_ends(row):
        try:
            parent_dates = mdf[mdf.autolin == row.parent].date
            child_dates = mdf[mdf.autolin == row.proposed_sublineage].date
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
    def get_regions(lin):
        try:
            # return ",".join(list(mdf.loc[[l for l in t.get_leaves_ids(nid) if l in mdf.index]].country.value_counts().index))
            return ",".join(list(mdf[mdf.autolin == lin].country.value_counts().index))
        except KeyError:
            return np.nan
    def get_regions_percents(lin):        
        try:
            # return ",".join([str(round(p,2)) for p in mdf.loc[[l for l in t.get_leaves_ids(nid) if l in mdf.index]].country.value_counts(normalize=True)])
            return ",".join([str(round(p,2)) for p in mdf[mdf.autolin == lin].country.value_counts(normalize=True)])
        except KeyError:
            return np.nan
    pdf['child_regions'] = pdf.proposed_sublineage.apply(get_regions)
    pdf['child_regions_count'] = pdf.child_regions.apply(lambda x:x.count(",")+1)
    # pdf['ParentRegions'] = pdf.parent_nid.apply(get_regions)
    pdf['child_region_percents'] = pdf.proposed_sublineage.apply(get_regions_percents)
    # pdf['ParentRegionPercents'] = pdf.parent_nid.apply(get_regions_percents)
    def host_jump(lin):
        try:
            return mdf[mdf.autolin == lin].host.nunique() > 1
            # return mdf.loc[[l for l in t.get_leaves_ids(proposed_sublineage_nid) if l in mdf.index]].host.nunique() > 1
        except:
            return False
    print("Identifying host jumps.")
    pdf['host_jump'] = pdf.proposed_sublineage.apply(host_jump,axis=1)
    print("Generating cov-spectrum URLs.")
    def generate_url(row):
        child_mset = t.get_haplotype(row.proposed_sublineage_nid)
        parent_mset = t.get_haplotype(row.parent_nid)
        net_mset = child_mset - parent_mset
        mset_str = "[{}-of:{}]".format(len(net_mset), ", ".join([m[1:] for m in net_mset]))
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
            td = (row.latest_child - row.earliest_child)
            time = (td.days-td.days%7)/7 + 1
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
                #further filter aa changes in orf1a/b so that they're properly processed for taxonium viewing and not counted redundantly
                #in our code, ORF1a changes are annotated as both ORF1a and ORF1ab, ORF1b are annotated as ORF1ab only.
                aas = []
                for a in translation.get(n.id,[]):
                    if a.original_aa == a.alternative_aa:
                        continue #ignore synonymous mutations
                    aas.append(a)
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
        pdf['taxlink'] = pdf.apply(lambda row:write_taxonium_url(row.parent, changes_to_list(row.mutations)),axis=1)
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