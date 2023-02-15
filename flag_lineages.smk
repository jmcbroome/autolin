import sys
from pathlib import Path
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

rule open_pull_request:
    input:
        "{tree}.proposed.report.tsv",
        "{tree}.proposed.pb",
        "{tree}_metadata.viz.tsv"
    output:
        "{tree}.pullreq.report.tsv"
    run:
        commandstr = "{config[python]} open_pull_request.py -r {config[request_params][designation_repo]} \
            -i {input[0]} -t {input[1]} -s {config[request_params][valid_samples]} -c {config[request_params][representative_number]} \
            -u {config[request_params][countries]} -m {config[request_params][maximum]} \
            -o {output} -a {config[request_params][active_since]} -d {config[request_params][samples_different]}"
        if eval(str(config["request_params"]["local_only"])):
            commandstr += " --local"
        if eval(str(config["request_params"]["auto_merge"])):
            commandstr += " --automerge"
        if eval(str(config["request_params"]["no_prefix"])):
            commandstr += " --no-prefix"
        if eval(str(config["request_params"]["growth_model"]["use_model"])):
            commandstr += " --model_growth \
                           --metadata {input[2]} \
                           --draws {config[request_params][growth_model][draws]} \
                           --tune {config[request_params][growth_model][tune]} \
                           --target_accept {config[request_params][growth_model][target_accept]} \
                           --min_country_weeks {config[request_params][growth_model][min_country_weeks]} \
                           --maxperc {config[request_params][growth_model][max_proportion_considered]}"
        shell(commandstr)

rule write_issues:
    input:
        "{tree}.proposed.report.tsv",
        "{tree}.proposed.pb",
        "{tree}_metadata.viz.tsv"
    output:
        "{tree}.issues.log"
    run:
        if eval(str(config["issue_params"]["local_only"])):
            shell("{config[python]} write_issues.py --local -i {input[0]} -t {input[1]} -m {input[2]} -n {config[issue_params][number]} -s {config[issue_params][sort_by]} -p {config[issue_params][prefix]} -c {config[issue_params][samples_named]}")
        else:
            shell("{config[python]} write_issues.py -i {input[0]} -t {input[1]} -m {input[2]} -n {config[issue_params][number]} -s {config[issue_params][sort_by]} -p {config[issue_params][prefix]} -c {config[issue_params][samples_named]}")

rule add_metadata:
    input:
        "{tree}.metadata.tsv",
        "{tree}.proposed.pb"
    output:
        "{tree}_metadata.viz.tsv"
    shell:
        "{config[python]} extract_annotations_viz.py {input[1]} {input[0]} {output}"

rule generate_report:
    input:
        "{tree}.proposed.pb",
        "{tree}.proposed.tsv",
        "{tree}.metadata.tsv"
    output:
        "{tree}.proposed.report.tsv"
    shell:
        "{config[python]} generate_lineage_report.py -i {input[0]} -p {input[1]} -o {output} -f {config[reference_genome]} -g {config[reference_gtf]} -m {input[2]} -r {config[lineage_params][downloadable_samples]}" #-d {config[lineage_params][earliest_date]}

rule propose:
    input:
        "{tree}.filtered.pb",
        "escape_weights.tsv",
        "{tree}.sample_weights.tsv"
    output:
        "{tree}.proposed.tsv",
        "{tree}.proposed.pb"
    log:
        "{tree}.proposal.log"
    run:
        d = {"input":input[0]}
        for k,v in config['lineage_params'].items():
            if k == 'weight_params' or k == 'earliest_date' or k == 'downloadable_samples':
                continue #these are used elsewhere.
            if v in ['None','True','False']:
                d[k] = eval(v)
            else:
                d[k] = v
        d['aaweights'] = input[1]
        d['samples'] = input[2]
        d['verbose'] = True
        d['reference'] = config['reference_genome']
        d['gtf'] = config['reference_gtf']
        d['dump'] = output[0]
        d['output'] = output[1]
        command_args = []
        for k,v in d.items():
            if type(v) == bool:
                if v:
                    command_args.append("--" + str(k))
            else:
                command_args.append("--" + str(k) + " " + str(v))
        print("FULL PROPOSAL COMMAND: ", "{config[python]} propose_sublineages.py " + " ".join(command_args))
        shell("{config[python]} propose_sublineages.py " + " ".join(command_args) + "> " + log[0])

rule unzip_metadata:
    input:
        "{tree}.metadata.tsv.gz"
    output:
        "{tree}.metadata.tsv"
    shell:
        "gunzip -c {input} > {output}"

rule compute_region_weights_and_dates:
    input:
        "{tree}.metadata.tsv"
    output:
        temp("{tree}.sample_weights.tsv")
    run:
        mdf = pd.read_csv(input[0],sep='\t')
        def get_dt(dstr):
            try:
                if type(dstr) == str:
                    return dt.datetime.strptime(dstr,"%Y-%m-%d")
                else:
                    return dstr
            except:
                return np.nan
        mdf['Date'] = mdf.date.apply(get_dt)
        mindate = get_dt(config['lineage_params']['earliest_date'])
        target = mdf[mdf.Date.dt.date >= mindate]
        scale = config['lineage_params']['weight_params']['country_weighting']
        invweights = 1/target.country.value_counts(normalize=True)
        to_use = (invweights-invweights.min())/(invweights.max()-invweights.min()) * scale + 1
        with open(output[0],"w+") as outf:
            for i,d in target.iterrows():
                print(d.strain + "\t" + str(to_use[d.country]),file=outf)        

rule compute_escape_weights:
    output:
        temp("escape_weights.tsv")
    run:
        with open("escape_weights.tsv","w+") as wout:
            edf = pd.read_csv(config['escape_data'])
            wvc = edf.groupby("label_site").site_total_escape.mean()
            for ls in wvc.index:
                print("\t".join([str(v) for v in ['S', ls[1:]+ls[0], 1 + wvc[ls] * config['lineage_params']['weight_params']['escape_weighting']]]), file=wout)

rule download_tree:
    #many public trees can be obtained online, from http://hgdownload.soe.ucsc.edu/goldenPath/wuhCor1/UShER_SARS-CoV-2/2022/09/28/
    #if the indicated tree file is not found locally, it will attempt to download it from this source.
    output:
        "{tree}.pb.gz",
        "{tree}.metadata.tsv.gz"
    run:
        date = wildcards.tree.split(".")[0].split("-")[1:]
        link = "http://hgdownload.soe.ucsc.edu/goldenPath/wuhCor1/UShER_SARS-CoV-2/{}/{}/{}/".format(*date)
        fn = wildcards.tree.split(".")[0] + ".all.masked.pb.gz"
        shell("wget " + link + fn)
        shell("mv {} {}".format(fn, output[0]))
        shell("wget " + link + output[1])

rule filter_tree:
    input:
        "{tree}.pb.gz",
        "{tree}.nodestats.txt"
    output:
        temp("{tree}.filtered.allan.pb")
    shell:
        "{config[python]} filter_reversions.py {input[0]} {input[1]} {output}"

rule simplify_annotations:
    input:
        "{tree}.filtered.allan.pb"
    output:
        temp("{tree}.filtered.pb")
    shell:
        "{config[python]} strip_annotations.py {input} {output}"

rule collect_node_statistics:
    input:
        "{tree}.pb.gz"
    output:
        "{tree}.nodestats.txt"
    conda:
        #usher is being handled in a separate conda environment due to ongoing environment conflicts over boost-cpp versions between UShER and earlier versions of BTE
        #at some point these should be resolved and UShER/matUtils can be added to the env.yml. 
        "usher.yml"
    shell:
        "matUtils summary -i {input} -N {output}"