# Aphids chemosensory genes comparative genomic

## Selection signatures in chemosensory genes

For each OR and GR families, CDS sequences were translated to protein sequences using Biopython https://biopython.org/docs/1.75/api/Bio.Seq.html, and aligned using MUSCLE https://drive5.com/muscle/. Phylogenetic tree was computed using protein sequences alignment and IQTREE http://www.iqtree.org/. 
Codon alignment were obtained from proteins alignment and CDS sequences using PAL2NAL perl script, with paml as output format.
Gaps were kept in codon alignments.

These codon alignments were used as input for codeml paml, with the branch and the site model. The branch model assumes different ω ratios for different branches on the phylogeny and was used to detect positive selection acting on particular sequences of the subfamilies without averaging ω throughout the phylogenetic tree. We also tested the site model that can be used to identify positively selected sites in a multiple-sequence alignment. 


## Aphids chemosensory genes synteny

To explore chemosensory genes synteny, we used R package syntenet https://github.com/almeidasilvaf/syntenet.
This analysis requires the gene annotation gff files and the files containing the protein sequences of the annotated genes, for all genes.
We also have genes of interest lists (here the chemosensory genes). 
To explore synteny, we need the tabular output of diamond https://github.com/bbuchfink/diamond

We also explored mode of duplication of chemosensory genes using R package doubletrouble https://github.com/almeidasilvaf/doubletrouble
which needs also tabular output of diamond. 

## TE enrichment around chemosensory genes 

To determine whether some transposable elements (TE) are found enriched around chemosensory genes, we used R package LOLA https://github.com/nsheff/LOLA. 

This analysis requires the gene annotation gff files for all genes and for genes of interest and the TE annotation gff files.
The genes gff were converted to bed format using AGAT agat_convert_sp_gff2bed.pl function. 
The genes coordinates were enlarged (2 kb and 10 kb) using bedtools slop https://bedtools.readthedocs.io/en/latest/content/tools/slop.html. 
Then the genes coordinates were sorted and merged. 

The TE annotation files were converted to the appropriate bed format using an homemade perl script. 

Once enriched TE list was obtained, we determined relationships between these TE and the genes of interest using an adapted TEgrip algorithm https://github.com/marieBvr/TEs_genes_relationship_pipeline. We also calculated divergence between TE copies, using Kimura distance and compared TE copies associated to chemosensory genes and all TE copies. 




