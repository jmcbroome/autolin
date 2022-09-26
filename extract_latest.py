import pandas as pd
import datetime as dt
import numpy as np
mdf = pd.read_csv("public-2022-08-29.metadata.tsv",sep='\t')
def get_dt(dstr):
    try:
        return dt.datetime.strptime(dstr,"%Y-%m-%d")
    except:
        return np.nan

mdf['Date'] = mdf.date.apply(get_dt)
rdf = pd.read_csv("late.proposed.lineages.report.tsv",sep='\t')
target = mdf[mdf.Date > dt.datetime(month=6,day=30,year=2022)]
scale = 100
invweights = 1/target.country.value_counts(normalize=True)
to_use = (invweights-invweights.min())/(invweights.max()-invweights.min()) * scale + 1
with open("later_samples.weights.txt","w+") as outf:
    for i,d in target.iterrows():
        print(d.strain + "\t" + str(invweights[d.country]),file=outf)