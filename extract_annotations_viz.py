import bte
import sys
import pandas as pd
t = bte.MATree(sys.argv[1])
df = pd.read_csv(sys.argv[2],sep='\t')
df['auto_annotation'] = df.strain.apply(lambda x:t.get_node(x).most_recent_annotation()[0])
df.to_csv(sys.argv[3],sep='\t',index=False)
