# automate-lineages-prototype

This repository contains an implementation of a concept for automated lineage annotation for parsimony phylogenetics, with a specific
emphasis on SARS-CoV-2. 

Robust viral lineage identification is a critical problem. Many epidemiological analyses are reliant on abstracting a viral phylogeny
into genotype categories, or lineages, such as Omicron and Delta. However, lineages are often identified in a haphazard fashion, with individual researchers scrutinizing the phylogeny manually and proposing new lineages individually, in a way informed by their own biases and focuses. Here, I propose a simple alternative lineage identification concept based on genotype representation.

Intuitively, to generate a single lineage annotation, we identify
the node on the tree with the highest scaled mean genotype representation. Stating that a sample is descended from a node conveys something about the genotype of that samples. For each sample descended from a given node, that node
represents some proportion of its genotype- that proportion being the number of mutations accumulated by the time of that ancestor as a 
proportion of the total number of mutations that sample eventually has. A given node often has many descendents, and it represents some
part of the genotype of each of them. A good candidate for a lineage label is a node that represents substantial genetic information across many descendents. We can therefore partition the tree with a new lineage at the node with the highest scaled mean genotype representation.

We define the scaled mean representation for a given internal node as 

$$
N * D \over{{S \over{N}} + D}
$$

where N is the number of descendent samples, S is the sum of all path distances (e.g. branch length separating a specific sample and the given internal node) to descendent samples, and D is the path distance from the given node to the tree root. Distance, in the case of a Mutation Annotated Tree, is in number of mutations accumulated. This value increases both as the number of descendents increase and as the relative distance from the root versus the distance to the descendents increases. 

However, choosing the maximum value of this metric only identifies a single node on the tree as a putative lineage. Once a single best lineage node has been identified, additional lineages are generated in two general ways. First, additional lineages 
can be generated serially, by ignoring all samples descended from a currently identified lineage and computing the maximum weight among 
nodes from the remaining samples, until most samples are included in a lineage or a fixed number have been generated. Secondly, additional
lineages can be generated hierarchically (as sublineages), which is performed by repeating the serial process on a subtree defined by 
a higher level lineage label. 

Critically, this system allows for flexible alterations of sample or mutational weight when computing these metrics. Mutations which are known to be associated with changes in epidemiological behavior can be given high weights, meaning they contribute a larger value to the path length and will push the criterion to define lineages that include them on their path. Similarly, samples from regions of interest or underrepresented times and places can be counted as more than one sample, again pushing the criterion to define lineages which represent these specific samples well.

This does not invalidate current approaches or lineages defined under other systems, but instead augments those lineage systems
by proposing new sublineages automatically. 

## Implementation

The prototype implementation in "propose_sublineages.py" is dependent on [BTE](https://github.com/jmcbroome/BTE). You can install BTE from bioconda.

```
conda install -c bioconda bte
```

This prototype is in very early stages and is likely buggy. It takes a tree in the MAT format, usually downloaded from [here](http://hgdownload.soe.ucsc.edu/goldenPath/wuhCor1/UShER_SARS-CoV-2/). This repository also includes a small test MAT containing only 
10,000 samples.

Currently available parameters are
- clear: remove all current annotations from the tree. If not used, new sublineages will be proposed where possible from current lineages.
- recursive: add additional sublineages to any new proposed sublineages where possible. If not used, sublineages to existing lineages will be generated only.
- dump: write a table containing proposed sublineages.
- output: write a MAT protobuf containing all proposed sublineages.
- labels: write a table containing samples and their outermost lineage assignments. Used for visualization with matUtils extract -j.
- floor: set the minimum mean representation gain to this value. Default 0
- maxpath: set the maximum Rd value, which is currently computed as the distance from the parent to the nearest grandparent lineage or root if no grandparent lineage exists.

## Current Issues

### Tree Growth and Lineage Optimization

This method is designed with a fixed tree, and designates lineages which are optimal representations given an existing tree. In reality, however, new data comes in and the tree grows; a formerly optimal lineage may no longer be optimal after the tree grows or is optimized itself. At the moment we have no good solution to the issue of the accumulation of suboptimal lineages as the tree grows, and are interested in feedback or suggestions for the effective identification of when a lineage should be retracted and reassigned. 

### Regularization of Lineage System Size

We're also interested in a method to effectively penalize the addition of new lineages. This system is an automatic way to identify the best choice for a new lineage label, but is not effective at deciding whether a new lineage label should be added at all. We currently use a variety of thresholds for minimum lineage size and distinction from its parent lineage, but ideally, we would have a more informed regularization scheme.

## Feedback

Feel free to try this method out! Please report bugs and give feedback in the issues. You can also email me at jmcbroom@ucsc.edu.
