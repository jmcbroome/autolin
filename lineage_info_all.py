import subprocess
import datetime as dt
import pandas as pd
from os.path import exists
end_date = dt.datetime(year=2022,month=8,day=24) #the latest date collected.
start_date = dt.datetime(year=2021,month=5,day=25) #the earliest date collected.
#proceed in reverse order.
# cdate = start_date
pdf = {k:[] for k in ['Lineage','Parent','DateAdded','EarliestDate','LatestDate','Count']}
dfv = []
for date in pd.date_range(start_date,end_date):
    datestr = date.strftime("%Y-%m-%d")
    #use a separate process so the memory use doesn't explode from loading many consecutive trees
    #check to see if this file already exists for convenience.
    if exists(datestr + ".linfo.tsv"):
        print("{} already exists; continuing".format(datestr + ".linfo.tsv"))
        continue
    try:
        subprocess.check_call("python3 automate-lineages-prototype/collect_lineage_info.py " + datestr,shell=True)
    except KeyboardInterrupt:
        print("User requested cancel. Exiting")
        exit(1)
    except:
        print("Something went wrong with day {}; continue".format(datestr))
        continue
    subdf = pd.read_csv(datestr+".linfo.tsv",sep='\t')
    subdf['DateCollected'] = datestr
    dfv.append(subdf)
pdf = pd.concat(dfv)
pdf.to_csv("linfo_merged.tsv",sep='\t',index=False)
