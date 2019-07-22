# EXACT
Exact inference under the perfect phylogeny model
Motivation: Many inference tools use the
perfect phylogeny model
(PPM) to learn trees from noisy
variant allele frequency
(VAF) data. Learning in this setting is hard. Existing tools are based on
approximate or heuristic algorithms, leaving learning under the PPM with noisy data still a poorly
understood task.
Results: We close this knowledge-gap in the scenario where the mutations that are relevant for
evolution can be clustered into a small number of groups, and the trees to be reconstructed have
a small number of nodes. We use a careful combination of algorithms, software, and hardware, to
develop EXACT: a tool that can explore the space of all possible phylogenetic trees, and perform
exact inference under the PPM with noisy data. EXACT allows users to obtain not just the most-
likely tree for some input data, but exact statistics about the distribution of trees that might explain
the data. We show that EXACT outperforms several existing tools for this same task.
