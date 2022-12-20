import pandas as pd
import numpy as np
import pymc3 as pm

def get_growth_model(df, min_data = 5, target_accept = 0.8, tune = 1000, draws = 1000):
    #prepare the dataset overall.
    rc = df.groupby(['country','auto_annotation',pd.Grouper(key='date', freq='1W')]).strain.count().reset_index()
    rc = rc.rename({"strain":"count"},axis=1).sort_values("date")
    cc = rc.groupby(["date","country"])['count'].sum().to_dict()
    rc['country_count'] = rc.apply(lambda row:cc.get((row.date,row.country)),axis=1)
    rc['country_perc'] = rc['count'] / rc.country_count
    #split the data by annotation.
    growd = {}
    for ann, osdf in rc.groupby("auto_annotation"):
        #prep the data.
        X_week = []
        X_country_week_total = []
        Y = []
        #start with one annotation as a test.
        for keys, sdf in osdf.groupby(['country']):
            sdf = sdf.sort_values('date').reset_index()
            #skip country/lineage pairs which don't have at least two weeks. Single weeks can only be used for initial proportion inference and are not generally useful.
            if sdf.shape[0] >= 2:
                for i,d in sdf.iterrows():
                    X_week.append(i)
                    X_country_week_total.append(d.country_count)
                    Y.append(d['count'])
        if len(Y) < min_data:
            print(f"Skipping {ann} for insufficient data ({len(Y)} consecutive countryweeks).")
            continue
        print(f"Fitting model on annotation {ann} with {len(Y)} consecutive countryweeks of data.")
        #fit the model!
        with pm.Model() as model:
            growth = pm.Normal(name='growth', sd=5)
            initial_proportion = pm.TruncatedNormal(name='initial_proportion',upper=1,lower=0)
            #cap the initial proportion value at 1 (100%) and log it for use. 
            #This value will be informed for week 0 of a set of values. It should vary across countries, but not by too much.
            log_initial_proportion = pm.Deterministic(name="log_initial_proportion",var=np.log(initial_proportion))
            #estimate our expected proportion for this week, given our initial proportion and week. correct it back by exponentiation.
            current_proportion = pm.Deterministic(name='base_proportion', var = np.e**(log_initial_proportion + growth * X_week))
            #sampling process with our actual observed values.
            y_obs = pm.Binomial(name='sampled', n=X_country_week_total, p=current_proportion, observed=Y)
            #perform the actual inference process.
            idata = pm.sample(draws=draws,tune=tune,progressbar=False,return_inferencedata=True,target_accept=target_accept)
        growd[ann] = np.percentile(idata.get_values("growth"), [5,50,95])
    return growd