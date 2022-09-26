import pandas as pd
edf = pd.read_csv("escape_data.csv")
wvc = edf.groupby("label_site").site_mean_escape.mean()
for ls in wvc.index:
    print("\t".join([str(v) for v in ['S', ls[1:]+ls[0], 1 + wvc[ls] * 100]]))