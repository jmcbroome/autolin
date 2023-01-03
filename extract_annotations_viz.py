import sys
sys.path.append("~/bin:")
import bte
import pandas as pd
t = bte.MATree(sys.argv[1])
df = pd.read_csv(sys.argv[2],sep='\t')
def pull_annotations(sample):
    try:
        n = t.get_node(sample)
    except ValueError:
        print("WARNING: Can't load node object and obtain annotations for {}".format(sample))
        return "None"
    return n.most_recent_annotation()[-1]
df['auto_annotation'] = df.strain.apply(pull_annotations)
df.to_csv(sys.argv[3],sep='\t',index=False)
