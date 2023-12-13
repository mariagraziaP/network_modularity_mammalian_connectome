# network_modularity_mammalian_connectome

## Details

Here you'll find MATLAB scripts and functions to analyze the mammalian connectomes from the MAMI database. Results are reported in the paper "Relation of connectome topology to brain volume across 103 mammalian species".

## Paper abstract

The brain connectome is an embedded network of anatomically interconnected brain regions and the study of its topological organization in mammals has become of paramount importance, due to its role in scaffolding brain function and behavior. Unlike many other observable networks, brain connections incur material and energetic cost, and their length and density are volumetrically constrained by the skull. Thus, an open question is how differences in brain volume impact connectome topology. We address this issue using the MaMI database, a diverse set of mammalian connectomes reconstructed from 201 animals, covering 103 species and 12 taxonomy orders, whose brain size varies over more than four orders of magnitude. Our analyses focus on relationships between volume and modular organization. After having identified modules through a multi-resolution approach, we observed how connectivity features relate to the modular structure and how these relations vary across brain volume. We found that as the brain volume increases modules become more spatially compact and dense, comprising more costly connections. Furthermore, we investigated how spatial embedding shapes network communication, finding that as brain volume increases, nodes' distance progressively impacts communication efficiency. We identified modes of variation in network communication policies, as smaller and bigger brains show higher efficiency in routing- and diffusion-based signaling, respectively. Finally, bridging network modularity and communication, we found that in larger brains modular structure imposes stronger constraints on network signaling. Altogether our results show that brain volume is systematically related to mammalian connectome topology and spatial embedding imposes tighter restrictions on larger brains.

## Libraries

This project relies on external libraries:
- BCT toolbox: to extract measures on brain networks (https://sites.google.com/site/bctnet/).
- georand: package to compute surrogate networks ("The contribution of geometry to the human connectome", J.A. Roberts et al., NeuroImage, 2016, DOI: https://doi.org/10.1016/j.neuroimage.2015.09.009).
