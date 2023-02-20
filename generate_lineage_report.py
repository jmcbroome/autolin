import sys
sys.path.append("./SARS2_RBD_Ab_escape_maps/")
import bindingcalculator as bc
import bte
import pandas as pd
import numpy as np
import datetime as dt
import argparse
from urllib import parse
import requests

def argparser():
    parser = argparse.ArgumentParser(description="Compute detailed lineage reports for all existing lineages in the tree.")
    parser.add_argument("-i", "--input", required=True, help='Path to protobuf to compute reports from.')
    parser.add_argument("-p", "--proposed", required=True, help='Path to the file containing dumped sublineage proposals.')
    parser.add_argument("-o", "--output", help='Name of the output table.',default=None,required=True)
    parser.add_argument("-m", "--metadata", help="Path to a metadata file matching the protobuf.",required=True)
    parser.add_argument("-f", "--reference", default=None, help="Path to a reference fasta file. Use with -g to annotate amino acid changes and immune escape in the expanded output.")
    parser.add_argument("-g", "--gtf", default=None, help="Path to a reference gtf file. Use with -f to annotate amino acid changes and immune escape in the expanded output.")
    parser.add_argument("-d", "--date", default=None, help="Ignore individual samples from before this date when computing reports. Format as %Y-%m-%d")
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
    final = urlbase + queries
    if len(final) > 2000:
        return "Taxonium link could not be generated- too many characters."
    return final

def update_aa_haplotype(caas, naas):
    #add all amino acid mutations in naas (new amino acids) to caas (current amino acids) if they don't overlap with an existing protein index.
    cindexes = set([(aa.gene, aa.aa_index) for aa in caas])
    for naa in naas:
        #ignore synonymous mutations
        if naa.original_aa != naa.alternative_aa:
            #we proceed out to in, so sites are "locked in" once a change is seen.
            if (naa.gene, naa.aa_index) not in cindexes:
                caas.append(naa)
    return caas

def fill_output_table(t,pdf,mdf,fa_file=None,gtf_file=None,mdate=None):
    print("Filling out metadata with terminal lineages.")
    def get_latest_lineage(s):
        for anc in t.rsearch(s):
            if len(anc.annotations) > 0:
                if len(anc.annotations[0]) > 0:
                    return anc.annotations[0]
    mdf['date'] = mdf.date.apply(get_date)
    if mdate != None:
        mdf = mdf[mdf.date > dt.datetime.strptime(mdate,"%Y-%m-%d")]
    mdf['autolin'] = mdf.strain.apply(get_latest_lineage)
    mdf.set_index('strain',inplace=True)
    #parent lineage size has to be inclusive to get a sensible percentage.
    def parent_lineage_size(lin):
        return mdf[mdf.pango_lineage_usher == lin].shape[0]
    print("Computing sublineage percentages.")
    pdf['parent_lineage_size'] = pdf.parent.apply(parent_lineage_size)
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
    def get_start_ends(row):
        try:
            parent_dates = mdf[mdf.pango_lineage_usher == row.parent].date.dropna()
            child_dates = mdf[mdf.autolin == row.proposed_sublineage].date.dropna()
            return min(parent_dates),max(parent_dates),min(child_dates),max(child_dates)
        except KeyboardInterrupt:
            raise KeyboardInterrupt
        except ValueError:
            return np.nan,np.nan,np.nan,np.nan    
    print("Computing start and end dates.")
    applied_pdf = pdf.apply(lambda row: get_start_ends(row), axis='columns', result_type='expand')
    pdf = pd.concat([pdf, applied_pdf], axis='columns')
    pdf = pdf.rename({0:'earliest_parent',1:'latest_parent',2:'earliest_child',3:'latest_child'},axis=1)
    pdf['log_score'] = np.log10(pdf.proposed_sublineage_score)
    print("Tracking country composition.")
    def get_regions(lin):
        try:
            return ",".join(list(mdf[mdf.autolin == lin].country.value_counts().index))
        except KeyError:
            return np.nan
    def get_regions_percents(lin):        
        try:
            return ",".join([str(round(p,2)) for p in mdf[mdf.autolin == lin].country.value_counts(normalize=True)])
        except KeyError:
            return np.nan
    pdf['child_regions'] = pdf.proposed_sublineage.apply(get_regions)
    pdf['child_regions_count'] = pdf.child_regions.apply(lambda x:x.count(",")+1)
    pdf['child_region_percents'] = pdf.proposed_sublineage.apply(get_regions_percents)
    def host_jump(lin):
        try:
            return mdf[mdf.autolin == lin].host.nunique() > 1
        except:
            return False
    print("Identifying host jumps.")
    pdf['host_jump'] = pdf.proposed_sublineage.apply(host_jump)
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
            # if n.id == row.parent_nid:
                # break
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
    pdf['growth_index'] = pdf.apply(get_growth_score,axis=1)
    if gtf_file != None and fa_file != None:
        print("Performing translation and computing antibody binding scores.")
        translation = t.translate(fasta_file = fa_file, gtf_file = gtf_file)
        calculator = bc.BindingCalculator(csv_or_url='SARS2_RBD_Ab_escape_maps/processed_data/escape_calculator_data.csv')
        hstrs = []
        cev = []
        pev = []
        nev = []
        for i, row in pdf.iterrows():
            child_aas = []
            parent_aas = []
            all_aas = []
            past_parent = False
            for n in t.rsearch(row.proposed_sublineage_nid,True):
                #further filter aa changes in orf1a/b so that they're properly processed for taxonium viewing and not counted redundantly
                #in our code, ORF1a changes are annotated as both ORF1a and ORF1ab, ORF1b are annotated as ORF1ab only.
                alla = translation.get(n.id,[])
                if n.id == row.parent_nid:
                    past_parent = True
                if past_parent:
                    #only mutations at or behind the parent node contribute to its haplotype.
                    parent_aas = update_aa_haplotype(parent_aas, alla)
                else:
                    child_aas = update_aa_haplotype(child_aas, alla)
            hstr = ",".join([aa.aa_string() for aa in child_aas])
            cspikes = [a.aa_index for a in all_aas if a.gene == 'S' and a.aa_index in calculator.sites and a.original_aa != a.alternative_aa]
            pspikes = [a.aa_index for a in parent_aas if a.gene == 'S' and a.aa_index in calculator.sites and a.original_aa != a.alternative_aa]
            child_escape = calculator.binding_retained(cspikes)
            parent_escape = calculator.binding_retained(pspikes)
            net_escape_gain = parent_escape - child_escape
            hstrs.append(hstr)
            cev.append(child_escape)
            pev.append(parent_escape)
            nev.append(net_escape_gain)
        pdf['aav'] = hstrs
        pdf['sublineage_escape'] = cev
        pdf['parent_escape'] = pev
        pdf['net_escape_gain'] = nev
    def get_reversions(subnid):
        reversions = []
        allm = set()
        past_parent = False
        for n in t.rsearch(subnid,True):
            if n.id != subnid:
                if any([len(a) > 0 for a in n.annotations]):
                    past_parent = True
            if past_parent:
                for m in n.mutations:
                    opposite = m[-1] + m[1:-1] + m[0]
                    if opposite in allm:
                        #record any mutations between the sublineage and the parent that are opposite of a parent mutation.
                        reversions.append(opposite)
            else:
                for m in n.mutations:
                    allm.add(m)
        if len(reversions) > 0:
            return ",".join(reversions)
        else:
            return "No Reversions"
    pdf['reversions'] = pdf.proposed_sublineage_nid.apply(get_reversions)
    def get_mset(mutations):
        mhap = []
        locs = set()
        for mset in reversed(mutations.split(">")):
            for m in mset.split(','):
                if len(m) > 0:
                    location = int(m[1:-1])
                    if location not in locs:
                        locs.add(location)
                        mhap.append(m[1:])
        return ','.join(mhap)
    pdf['mset'] = pdf.mutations.apply(get_mset)
    #remove any entries that have no mutations with respect to the parent.
    trackrev = open("reversion_proposals_blocked.log","w+")
    def log_mset(row):
        if len(row.mset) > 0:
            return True
        else:
            print(f"Proposal {row.proposed_sublineage} child of {row.parent} blocked for having no unique mutations; branch nid {row.proposed_sublineage_nid}, reversions {row.reversions}",file=trackrev)
            return False
    pdf = pdf[pdf.apply(log_mset,axis=1)]
    trackrev.close()
    def get_representative_download(row):
        #query on parent lineage + mutations instead
        #and use requests to see how many are available.
        check_query = f"https://lapis.cov-spectrum.org/open/v1/sample/aggregated?pangoLineage={row.parent}&nucMutations={row.mset}"
        response = requests.get(check_query)
        if response.status_code != requests.codes.ok:
            print(f"WARNING: Lapis Error Status Code {response.status_code} for link {check_query}")
            return np.nan
        elif response.json()['data'][0]['count'] == 0:
            print(f"No samples available for lineage proposal {row.proposed_sublineage}")
            return np.nan
        else:
            #return the fasta download version of this link.
            return f"https://lapis.cov-spectrum.org/open/v1/sample/fasta?pangoLineage={row.parent}&nucMutations={row.mset}"
    pdf['seqlink'] = pdf.apply(get_representative_download,axis=1)
    def get_epi_isls(row):
        #open version for if we ever have problems with the queries
        #query = f"https://lapis.cov-spectrum.org/open/v1/sample/gisaid-epi-isl?pangoLineage={row.parent}&nucMutations={row.mset}"
        query = f"https://lapis.cov-spectrum.org/gisaid/v1/sample/gisaid-epi-isl?pangoLineage={row.parent}&nucMutations={row.mset}&accessKey=9Cb3CqmrFnVjO3XCxQLO6gUnKPd"
        return query
    pdf['epi_isls'] = pdf.apply(get_epi_isls,axis=1)
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
    odf = fill_output_table(t,pdf,mdf,args.reference,args.gtf,args.date)
    odf.to_csv(args.output,sep='\t',index=False)

if __name__ == "__main__":
    main()