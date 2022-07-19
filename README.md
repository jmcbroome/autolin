# automate-lineages-prototype

This repository contains an implementation of a concept for automated lineage annotation for parsimony phylogenetics, with a specific
emphasis on SARS-CoV-2. 

Robust viral lineage identification is a critical problem. Many epidemiological analyses are reliant on abstracting a viral phylogeny
into genotype categories, or lineages, such as Omicron and Delta. However, lineages are often identified in a haphazard fashion, with individual researchers scrutinizing the phylogeny manually and proposing new lineages individually, in a way informed by their own biases and focuses. Here, I propose a simple alternative lineage identification concept focused on genotype representation.

Intuitively, to generate a single lineage annotation, we identify
the node on the tree with the highest scaled mean genotype representation. Stating that a sample is descended from a node conveys something about the genotype of that samples. For each sample descended from a given node, that node
represents some proportion of its genotype- that proportion being the number of mutations accumulated by the time of that ancestor as a 
proportion of the total number of mutations that sample eventually has. A given node often has many descendents, and it represents some
part of the genotype of each of them. A good candidate for a lineage label is a node that represents substantial genetic information across many descendents. We can therefore partition the tree with a new lineage at the node with the highest scaled mean genotype representation.

We define the scaled mean representation for a given internal node as $N * Sd \over{Ud \over{N} + Sd}$ where N is the number of descendent samples, Ud is the total path distance to descendent samples, and Sd is the distance from the given node to the tree root. Distance, in the case of a Mutation Annotated Tree, is in number of mutations accumulated. This value increases both as the number of descendents increase and as the relative distance from the root versus the distance to the descendents increases. 

However, choosing the maximum value of this metric only identifies a single node on the tree as a putative lineage. Once a single best lineage node has been identified, additional lineages are generated in two general ways. First, additional lineages 
can be generated serially, by ignoring all samples descended from a currently identified lineage and computing the maximum weight among 
nodes from the remaining samples, until most samples are included in a lineage or a fixed number have been generated. Secondly, additional
lineages can be generated hierarchically (as sublineages), which is performed by repeating the serial process on a subtree defined by 
a higher level lineage label. 

By itself, this index does not have a stopping condition or a way to evaluate the quality of a proposed sublineage in general. To address this critical weakness, we define the following:

$$
Sd \over{Ud \over{N} + Sd} - Rd \over{{Ud \over{N} + Sd + Rd} > M
$$

This equation evaluates whether a proposed sublineages representation (left) is greater than the theoretical representation of the root/parent lineage against a hypothetical grandparent lineage with Rd (or more) mutations separating them. Generally, these will be large positive values when the sublineage is a much better representative of its samples than the parent lineage, and negative when new sublineage labels represent less than Rd mutations per descendent sample compared to their parent. 

We can increase the minimum value (M) of the equation from the default of 0 and the Rd value from the default of 1. Both of these adjustments will prevent the creation of small, marginal sublineages that are barely better representations than existing labels.

With this metric, lineages are generated serially and while disregarding all samples labeled at this level until all remaining candidate nodes fail the inequality. Lineages are then generated hierarchically by individually generating sublineages at each outermost lineage (a lineage with no sublineages) until no candidate nodes pass the inequality. Sublineage proposal is complete when all existing outermost lineages have no descendent nodes with positive representation gains under the given parameters.

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

### Non-Tree Considerations

Many lineage systems take into account geographical or temporal elements. This is because many lineage systems do not strictly
seek to optimize genotype representation, but rather epidemiological relevance- while genotype is clearly critical to epidemiological behavior of a virus, other elements like location and relative sampling rate can affect this. We therefore want to optionally support
the use of metadata frequencies or other considerations in sample contribution to optimal lineage defintions. For example, we may want to more heavily weight samples from an underrepresented and undersampled area for lineage calls.

For example, we could weight samples by the inverse of their representation across the tree (if Peru consists of 5% of samples, then we can give Peruvian samples a weight multiplier
1/0.05 = 20). 

We might also consider some categories of mutation to be substantially more important than others. Critical spike protein mutations can lead to large epidemiological changes while being a relatively small part of a genotype. We can treat these mutations as contributing more than 1 to the total path length as we compute our genotype representation.

There are many potential systems, but which systems to consider directly supporting in any final implementation of this method is a critical question, though the basic method should remain based on phylogenetic information alone.

## Feedback

Feel free to try this method out! Please report bugs and give feedback in the issues. You can also email me at jmcbroom@ucsc.edu.
