# EXACT
Many inference tools use the perfect phylogeny model (PPM) to learn trees from noisy variant allele frequency (VAF) data. Learning in this setting is hard. Existing tools are based on approximate or heuristic algorithms, leaving learning under the PPM with noisy data still a poorly understood task. We close this knowledge-gap in the scenario where the mutations that are relevant for evolution can be clustered into a small number of groups, and the trees to be reconstructed have a small number of nodes. We use a careful combination of algorithms, software, and hardware, to develop EXACT: a tool that can explore the space of all possible phylogenetic trees, and perform exact inference under the PPM with noisy data. EXACT allows users to obtain not just the most-likely tree for some input data, but exact statistics about the distribution of trees that might explain the data. We show that EXACT outperforms several existing tools for this same task.

This repository includes code for EXACT, as well as code, and Matlab wrappers, for four competing tools: AncesTree, Canopy, CITUP, and PhyloWGS. There is a folder for each tool, include for EXACT, inside of which there is a run_me.m Matlab file that illustrates the use of each tool. This repository includes binaries for each tool that have been compiled for a Linux Ubuntu distribution and Intel IA64 architecture.

We include source code for tool inside each respective folder. The code for the four competing tools was obtained from their respective repositories, and are included here for convinience (https://github.com/raphael-group/AncesTreel; https://github.com/yuchaojiang/Canopy; https://github.com/sfu-compbio/citup; https://github.com/morrislab/phylowgs)



EXACT is under constant development. As a result, some features might be temporarily unavailable.
