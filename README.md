# automate-lineages-prototype

This repository contains an implementation of a concept for automated lineage annotation for parsimony phylogenetics, with a specific
emphasis on SARS-CoV-2. The concept is based on sum genotype representation- that is, an optimal lineage is one which represents as much
total genotype information across all of its descendents as possible. Ideally and in general, referring to a sample by its lineage should 
convey as much information about the genotype of that sample as possible- if I refer to a sample as an "Omicron" case, that 
conveys substantial information about the genotype of that sample. 

Robust viral lineage identification is a critical problem. Many epidemiological analyses are reliant on abstracting a viral phylogeny
into genotype categories, or lineages, such as Omicron and Delta. However, lineages are often identified in a haphazard fashion, with individual researchers scrutinizing the phylogeny manually and proposing new lineages individually, in a way informed by their own biases and focuses. Here, I propose a simple alternative lineage identification concept focused on genotype representation.

A set of slides overviewing the concept can be found [here](https://docs.google.com/presentation/d/1IrP0Ey2Ue1dtx2L-azNY1jLF4n-re7mB5BBM_m0327w/edit?usp=sharing). To summarize, to generate a single lineage annotation, we identify
the node on the tree with the highest sum genotype representation. Stating that a sample is descended from a node conveys something about the genotype of that samples. Intuitively, for each sample descended from a given node, that node
represents some proportion of its genotype- that proportion being the number of mutations accumulated by the time of that ancestor as a 
proportion of the total number of mutations that sample eventually has. A given node often has many descendents, and it represents some
part of the genotype of each of them. 

Once a single best lineage node has been identified, additional lineages are generated in two general ways. First, additional lineages 
can be generated serially, by ignoring all samples descended from a currently identified lineage and computing the maximum weight among 
nodes from the remaining samples, until most samples are included in a lineage or a fixed number have been generated. Secondly, additional
lineages can be generated hierarchically (as sublineages), which is performed by repeating the serial process on a subtree defined by 
a higher level lineage label. A fixed number of levels can be generated.

This does not invalidate current approaches or lineages defined under other systems, but can augment those lineage systems
by proposing new sublineages automatically.

## Current Issues

### Limiting Total Lineages Generated Intelligently

When to stop generating new lineages is unclear. If the goal is to maximize total representation of lineages and unlimited numbers can be identified, it will eventually degenerate to being exactly the phylogeny with one lineage per node and fail to be the phylogenetic summary we desire. We can cut off after a fixed number have been generated, but it would be preferable to have an optimization criteria that penalizes simply adding more and more lineage labels. I'd appreciate any suggestions along these lines.

### Online Phylogenetics

We also need to consider how this method could fit into an online framework, where we are adding to a preexisting lineage system.
We can extend the hierarchical approach and generate new lineages as sublineages, but deciding when a lineage needs to be subdivided
is less clear. We could simply set a maximum size using the size parameter and ensure no terminal lineage has more than X samples, but
surely there is a more informed approach that takes into account the representational power of a given lineage, and the potential gains
for further partitioning a preexisting lineage. This is closely related to the issue above.

### Non-Genetic Considerations

Many lineage systems take into account geographical or temporal elements. This is because many lineage systems do not strictly
seek to optimize genotype representation, but rather epidemiological relevance- while genotype is clearly critical to epidemiological behavior of a virus, other elements like location and relative sampling rate can affect this. We therefore want to optionally support
the use of metadata frequencies or other considerations in sample contribution to optimal lineage defintions. For example, we may want to more heavily weight samples from an underrepresented and undersampled area for lineage calls.

We currently provide a system for sample-level weighting, but providing prebuilt versions may be of value. For example, we could weight
samples by the inverse of their representation across the tree (if Peru consists of 5% of samples, then we can give Peruvian samples a weight multiplier
1/0.05 = 20). There are many potential systems, but which systems and what non-genetic information to consider directly supporting in any final implementation of this method is a critical question, though the basic method should remain based on genetic information alone.

## Implementation

The prototype implementation in "automate_lineages.py" is dependent on [BTE](https://github.com/jmcbroome/BTE). It is in very early stages and is likely buggy.

It takes a tree in the MAT format, usually downloaded from [here](http://hgdownload.soe.ucsc.edu/goldenPath/wuhCor1/UShER_SARS-CoV-2/). This repository also includes a small test MAT containing only 
10,000 samples.

Currently available parameters are
- count: generate at most this many annotations serially (mutually exclusive lineages)
- levels: generate this many levels of hierarchical annotations (inclusive lineages where one is a subset of another).
- size: the maximum number of samples to be included in a final level. If a lineage contains less than this number of samples at a level above the maximum indicated, that lineage will not be further subdivided.
- threshold: a minimum percentage of samples to cover with annotations at each level of annotation. A level is complete when either 
the maximum number of annotations are generated or the current set of annotations covers the threshold of samples.
- range: a number of index ranges to consider mutations within. Used to subset the genome to only consider more relevant regions,
e.g. the SARS-CoV-2 spike protein.
- weights: a sample-level multiplier to apply to node weights. For example, if we want to call lineages and specifically give weights 
for samples from the UK less weight, since they are heavily sequenced, we can apply a multiplier of 0.5 to the weights contributed by
samples from the UK to make them matter half as much as other samples on the tree for deciding optimal lineage locations.

## Feedback

Feel free to try this method out! You can install BTE from bioconda.

```
conda install -c bioconda bte
```

Please report bugs and give feedback in the issues. You can also email me at jmcbroom@ucsc.edu.