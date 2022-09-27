import sys
from pathlib import Path
sys.path.insert(0, Path(workflow.basedir).parent.as_posix())
import bte
from propose_sublineages import propose, argparser
import pandas as pd
import datetime as dt
import numpy as np

configfile: "config.yaml"

rule all:
    input:
        "{tree}.jsonl.gz"
    
rule taxonium:
    input:
        "{tree}.proposed.pb",
        "{tree}_metadata.viz.tsv"
    output:
        "{tree}.jsonl.gz"
    shell:
        "usher_to_taxonium -i {input[0]} -o {output} -m {input[1]} -c auto_annotation,pango_lineage_usher,country -g {config[genbank]}"

rule add_metadata:
    input:
        "{tree}.metadata.tsv",
        "{tree}.proposed.pb"
    output:
        "{tree}_metadata.viz.tsv"
    shell:
        #code that loads trees with BTE lives in individual scripts to ensure proper release of memory dedicated to large trees.
        "python3 extract_annotations_viz.py {input[1]} {input[0]} {output}"

rule generate_report:
    input:
        "{tree}.proposed.pb",
        "{tree}.proposed.tsv"
    output:
        "{tree}.proposed.report.tsv"
    shell:
        "python3 generate_lineage_report.py -i {input[0]} -p {input[1]} -o {output} -f {config[reference_genome]} -g {config[reference_gtf]}"
rule propose:
    input:
        "{tree}.pb.gz",
        "escape_weights.tsv",
        "{tree}.sample_weights.tsv"
    output:
        "{tree}.proposed.tsv",
        "{tree}.proposed.pb"
    log:
        "{tree}.log"
    run:
        args = argparser().parse_args(['--input',input[0]])
        d = vars(args)
        for k,v in config['lineage_params'].items():
            if k not in d.keys():
                continue
            elif v == "None":
                d[k] = None
            else:
                d[k] = v
        # d.update({k:v if v != "None" else k:None for k,v in  if k in d.keys()})
        d['aaweights'] = input[1]
        d['samples'] = input[2]
        d['verbose'] = True
        d['reference'] = config['reference_genome']
        d['gtf'] = config['reference_gtf']
        d['dump'] = output[0]
        d['output'] = output[1]
        propose(args)

rule unzip_metadata:
    input:
        "{tree}.metadata.tsv.gz"
    output:
        "{tree}.metadata.tsv"
    shell:
        "gunzip -c {input} > {output}"

rule compute_region_weights:
    input:
        "{tree}.metadata.tsv"
    output:
        "{tree}.sample_weights.tsv"
    run:
        mdf = pd.read_csv(input[0],sep='\t')
        def get_dt(dstr):
            try:
                return dt.datetime.strptime(dstr,"%Y-%m-%d")
            except:
                return np.nan
        mdf['Date'] = mdf.date.apply(get_dt)
        target = mdf[mdf.Date >= dt.datetime.strptime(config['lineage_params']['earliest_date'])]
        scale = config['lineage_params']['weight_params']['country_weighting']
        invweights = 1/target.country.value_counts(normalize=True)
        to_use = (invweights-invweights.min())/(invweights.max()-invweights.min()) * scale + 1
        with open(output[0],"w+") as outf:
            for i,d in target.iterrows():
                print(d.strain + "\t" + str(invweights[d.country]),file=outf)        

rule compute_escape_weights:
    output:
        "escape_weights.tsv"
    run:
        with open("escape_weights.tsv","w+") as wout:
            edf = pd.read_csv(config['escape_data'])
            wvc = edf.groupby("label_site").site_total_escape.mean()
            for ls in wvc.index:
                print("\t".join([str(v) for v in ['S', ls[1:]+ls[0], 1 + wvc[ls] * config['lineage_params']['weight_params']['escape_weighting']]]), file=wout)
