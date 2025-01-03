# BiocManager::install("almeidasilvaf/syntenet", force = T)
# BiocManager::install("almeidasilvaf/doubletrouble")
# BiocManager::install("ggtree")

set.seed(123) # for reproducibility
# Loading packages
library(syntenet)
library(doubletrouble)
library(here)
library(tidyverse)
library(ggtree)
library(patchwork)
library(readxl)
library(Biostrings)
library(reprex)

rm(list=ls())
here()
# dss2df <- function(dss) data.frame(width=width(dss), seq=as.character(dss), names=names(dss))

# fs::dir_tree(here("data"))


## GR and OR genes list 

# Iterate through sheets and read them as a list of tibbles
## GRs
gr_path <- here("data/data", "GR_list.xlsx")
gr_list <- gr_path %>%
    excel_sheets() %>%
    set_names() %>%
    map(read_excel, path = gr_path)

## ORs
or_path <- here("data/data", "OR_list.xlsx")
or_list <- or_path %>%
    excel_sheets() %>%
    set_names() %>%
    map(read_excel, path = or_path)


GR_list = do.call(rbind, lapply(seq_along(gr_list), function(i) {
    tmp = gr_list[[i]] %>% as.data.frame() 
})) %>% 
    mutate(Coordinate = gsub("\\(.*\\)", "", Coordinate)) %>% 
    mutate(Coordinate = gsub("\\.\\.", "-", Coordinate)) %>% 
    separate(Coordinate, into = c("contig", "coordinates"), sep = ":") %>% 
    separate(coordinates, into = c("start", "end"), sep = "-") %>% 
    mutate(start = as.numeric(start), end = as.numeric(end)) %>% 
    mutate(length = end - start + 1) %>% 
    mutate(species = substr(`Gene ID`,1,4)) %>% 
    mutate(geneSet = "GR")

GR_count = GR_list %>% 
    filter(! is.na(species)) %>% 
    group_by(species) %>% 
    tally()


OR_list = do.call(rbind, lapply(seq_along(or_list), function(i) {
    tmp = or_list[[i]] %>% as.data.frame()
})) %>% 
    mutate(Coordinate = gsub("\\(.*\\)", "", Coordinate)) %>% 
    mutate(Coordinate = gsub("\\.\\.", "-", Coordinate)) %>% 
    separate(Coordinate, into = c("contig", "coordinates"), sep = ":") %>% 
    separate(coordinates, into = c("start", "end"), sep = "-") %>% 
    mutate(start = as.numeric(start), end = as.numeric(end)) %>% 
    mutate(length = end - start + 1) %>% 
    mutate(species = substr(`Gene ID`,1,4)) %>% 
    mutate(geneSet = "OR")

liste_pseudo = c("Dpla_OR26","Dpla_OR32","Dpla_OR35","Dpla_OR2","Dpla_OR4","Dpla_OR7","Dpla_OR9","Dpla_OR10","Dpla_OR11","Dpla_OR12","Dpla_OR13","Dpla_OR14","Dpla_OR15","Dpla_OR16","Dpla_OR17","Dpla_OR18","Dpla_OR19","Dpla_OR20","Dpla_OR23","Dpla_OR24","Dpla_OR25","Dpla_OR31","Dpla_OR33","DplaGR9","DplaGR30","DplaGR31","DplaGR33","DplaGR40","Dvit_OR6","Dvit_OR15","Dvit_OR53","Dvit_OR13_B","Dvit_OR55","Dvit_OR58","Dvit_OR30","Dvit_OR9","Dvit_OR8","Dvit_OR11","Dvit_OR13_A","Dvit_OR26","Dvit_GR3","Dvit_GR18","Dvit_OR59","Dvit_OR57","Mcer_OR11","Mcer_OR20","Mcer_OR47","Rpad_OR11","Rpad_OR36","Rpad_OR52","Rpad_GR1","Rpad_OR10","Rpad_OR23","Rpad_OR49")

Chemo_list = rbind(OR_list, GR_list) %>% 
    mutate(pseudo = ifelse(`Gene ID` %in% liste_pseudo, "yes", "no")) 


OR_count = OR_list %>% 
    filter(! is.na(species)) %>% 
    group_by(species) %>% 
    tally()



## Load data

# Read GFF files as a GRangesList object
annot_files <- dir(here("data/data", "annotation"), full.names = TRUE)
annot <- lapply(annot_files, function(x) {
    granges <- rtracklayer::import(x, feature.type = "gene")
    ## Keep only useful fields
    granges <- granges[, c("source", "type", "ID", "Name")]
    return(granges)
})
names(annot) <- gsub(".gff", "", basename(annot_files))
annot <- GenomicRanges::GRangesList(annot)

# Read FASTA files as a list of AAStringSet objects
seq <- fasta2AAStringSetlist(here("data/data", "sequences"))

check_input(seq, annot)
# If check input is ok, then process data 
pdata <- process_input(seq, annot)


## You can have some differences on genes names in gff and in fasta, that have to be removed

# If check input is not ok, you can compare number of genes in seq and annot, they should be identical

ngenes <- data.frame(
    Genes = lengths(annot),
    Sequences = lengths(seq)
)

# Removed unwanted prefixes and suffixes from sequence names
seq <- lapply(seq, function(x) {
    names(x) <- gsub("-PA", "", names(x))
    names(x) <- gsub("-RA", "", names(x))
    names(x) <- gsub("-PB", "", names(x))
    names(x) <- gsub(".*gene=", "", names(x))
    names(x) <- gsub("\\.p\\d+", "", names(x))
    names(x) <- gsub("\\.t\\d+", "", names(x))
    return(x)
})

# Replace transcript IDs with gene IDs in sequence headers
# Get a list of data frames containing transcript-to-gene mapping

tx2gene <- lapply(annot_files, function(x) {
    granges <- rtracklayer::import(x)
    tgene <- as.data.frame(granges[granges$type == "mRNA"])
    if(nrow(tgene) < 1000) {
        tgene2 <- as.data.frame(granges[granges$type == "transcript"])
        tgene <- rbind(tgene, tgene2)
    }
    tgene_df <- data.frame(
        tx = as.character(tgene$ID),
        gene = as.character(tgene$Parent)
    )
    return(tgene_df)
})

names(tx2gene) <- gsub(".gff", "", basename(annot_files))

# Replace IDs
seq <- lapply(seq_along(seq), function(x) {
    tgene <- tx2gene[[x]]
    s1 <- seq[[x]][names(seq[[x]]) %in% tgene$tx]
    s2 <- seq[[x]][!names(seq[[x]]) %in% tgene$tx]
    names(s1) <- tgene[match(names(s1), tgene$tx), "gene"]
    s <- c(s1, s2)
    return(s)
})
names(seq) <- names(tx2gene)


# Replace protein IDs with gene IDs in sequence headers
# Get a list of data frames containing protein-to-gene mapping

tx2gene <- lapply(annot_files, function(x) {
    granges <- rtracklayer::import(x)
    tgene <- as.data.frame(granges[granges$type == "CDS"])
    tgene_df <- data.frame(
        tx = as.character(tgene$Name),
        gene = paste0("gene-", as.character(tgene$gene))
    )
    return(tgene_df)
})

names(tx2gene) <- gsub(".gff", "", basename(annot_files))

# Replace IDs

for (i in c(2,3,4,11)) {
    seq[i] <- lapply(seq[i], function(x) {
        tgene <- tx2gene[[i]]
        s1 <- seq[[i]][names(seq[[i]]) %in% tgene$tx]
        s2 <- seq[[i]][!names(seq[[i]]) %in% tgene$tx]
        names(s1) <- tgene[match(names(s1), tgene$tx), "gene"]
        s <- c(s1, s2)
        return(s)
    })
}

# Keep only longest sequence in case of duplication
seq <- lapply(seq, function(x) {
    s <- x[order(Biostrings::width(x), decreasing = TRUE), ]
    s <- s[!duplicated(names(s)), ]
    return(s)
})

# Rename column 'Name' in annotation to 'ID'
annot <- lapply(annot, function(x) {
    x$gene_id <- x$ID
    x$ID <- NULL
    return(x)
})




## Process data and analysis

# If diamond can not be run locally (like on windows), you have to export processed data to run diamond on another system.

export_sequences(pdata$seq, outdir = "./results/processed/all_genes")

# then you can read the diamond outputs 
diamond_list = read_diamond("./results/diamond_processed/all_genes")

# the function infer_syntenet detect synteny and infer the synteny network
net <- infer_syntenet(diamond_list, pdata$annotation, 
                      anchors = 5, max_gaps = 25, 
                      is_pairwise = F, verbose = T, 
                      outdir = "./results/infer_syntenet/all_genes")

# the function cluster_network allow to cluster the network in order to obtain the phylogenomic profiles 
clusters <- cluster_network(net)
profiles <- phylogenomic_profile(clusters)


# Create a vector with order of species to plot based on species tree
tree_order <- c(
    "MCE", "MPE", "DPL", "DNO", "API", "SMI", "PNI",
    "AGL", "AGO", "RMA", "RPA", "ELA", "DVI"
)

# Create a data frame of species metadata
species_metadata <- read_excel(here("data", "species_metadata.xlsx")) |>
    janitor::clean_names() |>
    filter(!is.na(acronyms)) |>
    mutate(acronyms = str_sub(acronyms, start = 1, end = 3)) |>
    dplyr::select(
        species, acronyms, origin, life_cycle, reproduction,
        host_specialization
    ) |>
    as.data.frame()
species_metadata <- inner_join(
    data.frame(acronyms = tree_order), species_metadata
)


# Add full names on names of order vector
species_order <- setNames(
    species_metadata$acronyms,
    species_metadata$species
)

# Exploring syntenic relationships among chemosensory genes

# Now, let’s filter our clusters data set to only contain clusters containing chemosensory genes.

# Combining OR and GR for each species into 2 character vectors
gr_genes <- unique(Reduce(rbind, gr_list)$`Gene ID`)
or_genes <- unique(Reduce(rbind, or_list)$`Gene ID`)
# Get clusters containing GR and OR genes
gr_clusters <- clusters |>
    mutate(Gene = str_replace_all(Gene, "[a-zA-Z]{3,5}_", "")) |>
    dplyr::filter(Gene %in% gr_genes)
or_clusters <- clusters |>
    mutate(Gene = str_replace_all(Gene, "[a-zA-Z]{3,5}_", "")) |>
    dplyr::filter(Gene %in% or_genes)

# Are there GR and OR in the same clusters (and thus, that are syntenic)

gr_or_clusters <- gr_clusters |>
    dplyr::rename(Gene_GR = Gene) |>
    dplyr::inner_join(or_clusters) |>
    dplyr::rename(Gene_OR = Gene)
nrow(gr_or_clusters)

#Interestingly, there are OR genes with syntenic relationships among themselves, and GR genes
#with syntenic relationships among themselves, but none of the GR genes are syntenic with
#OR genes. Synteny only occurs within the same class of chemosensory genes. Now,
#let’s see how many OR and GR genes are in synteny clusters.

# Get frequency of OR and GR genes that are syntenic
chemosensory_syntenic <- data.frame(
    GR_syntenic = nrow(gr_clusters),
    GR_syntenic_percentage = nrow(gr_clusters) / length(gr_genes),
    OR_syntenic = nrow(or_clusters),
    OR_syntenic_percentage = nrow(or_clusters) / length(or_genes)
)

#Now, let’s plot the phylogenomic profiles of these syntenic genes as a heatmap.

# Filter profile matrix to only keep clusters with syntenic OR and GR genes
profile_gr <- profiles[unique(gr_clusters$Cluster), ]
profile_or <- profiles[unique(or_clusters$Cluster), ]
# Plot profiles
## GR genes
p_prof_gr <- plot_profiles(
    profile_gr,
    # color by `host_specialization` and rename it to `specialization`
    species_annotation = species_metadata |>
        dplyr::select(acronyms, specialization = origin),
    cluster_species = species_order,
    show_colnames = TRUE,
    main = "Phylogenomic profiles of GR-encoding genes"
)

## OR genes
p_prof_or <- plot_profiles(
    profile_or,
    # color by `host_specialization` and rename it to `specialization`
    species_annotation = species_metadata  %>% 
        dplyr::select(acronyms, specialization = origin),
    cluster_species = species_order,
    show_colnames = TRUE,
    main = "Phylogenomic profiles of OR-encoding genes"
)


gs_clusters = find_GS_clusters(profile_gr, species_metadata %>% select(acronyms, host_specialization))
gs_clusters = find_GS_clusters(profile_or, species_metadata %>% select(acronyms, host_specialization))

id <- gs_clusters$Cluster[]
plot_network(net, clusters, cluster_id = id)
```

# Then wu use the doubletrouble R package to classify gene pairs

diamond_list = read_diamond("./results/diamond_processed/all_genes")

# Keep only intraspecies DIAMOND results
intra_comparisons <- paste0(
    names(pdata$annotation),
    "_",
    names(pdata$annotation)
)

diamond_intra <- diamond_list[intra_comparisons]
names(diamond_intra)
# Identifying and classifying duplicate pairs
classified_pairs <- classify_gene_pairs(
    blast_list = diamond_intra,
    annotation = pdata$annotation,
    scheme = "standard"
)


# Classifying genes

classified_genes <- classify_genes(classified_pairs)


toMatch = c("Or", "OR", "GR", "Gr")

Proximal_duplicated_genes = do.call(rbind, lapply(seq_along(classified_genes), function(i) {
    tmp = classified_genes[[i]] %>% as.data.frame() %>% 
        mutate(genome = names(classified_genes)[i]) %>% 
        filter(type == "PD") %>% 
        select(-type)
    tmp2 = tmp[grep(paste(toMatch,collapse="|"), x = tmp$gene, value = FALSE),]
    return(tmp2)
}))
write.table(Proximal_duplicated_genes, "./results/Proximal_duplicated_genes.csv", sep = "\t", row.names = FALSE)


# Create a species tree to include in the plot (for phylogenetic context)
b1 <- "(((((Mcerasi,Mpersicae),Dplantaginea),Dnoxia),(Apisum,Smiscanthi)),Pnigronervosa)"
b2 <- "((Aglycines,Agossypii),(Rmaidis,Rpadi))"
species_tree <- paste0(
    "(", b1, ",", b2, ")"
)
species_tree <- paste0(
    "(", species_tree, ",Elanigerum)"
)
species_tree <- paste0(
    "(", species_tree, ",Dvitifoliae);"
)
species_tree <- treeio::read.tree(text = species_tree)
# Plot species tree and get a vector of ordered species based on the phylogeny
p_tree <- ggtree(species_tree) +
    geom_tiplab() +
    theme_tree() +
    xlim(0, 17)
species_topology <- get_taxa_name(p_tree)

# Plot percentages of OG and GR genes by mode for each species
p_chemosensory_percentage <- Reduce(rbind, lapply(seq_along(classified_genes), function(x) {
    species <- names(classified_genes)[x]
    classified_genes[[x]]$species <- species
    return(classified_genes[[x]])
})) |>
    # Keep only OR and GR genes
    filter(!str_detect(gene, "gene")) |>
    mutate(
        class = case_when(
            str_detect(gene, "(?i)OR") ~ "OR",
            str_detect(gene, "(?i)GR") ~ "GR"
        )
    ) |>
    filter(!is.na(class)) |>
    dplyr::rename(Mode = type) |>
    group_by(species, class) |>
    dplyr::count(Mode) |>
    mutate(percentage = (n / sum(n) * 100)) |>
    ungroup() |>
    # Rename levels of column `type` to complete names
    mutate(
        Mode = str_replace_all(
            Mode,
            c(
                "DD" = "Dispersed",
                "PD" = "Proximal",
                "TD" = "Tandem",
                "WGD" = "Segmental"
            )
        )
    ) |>
    # Rename species
    mutate(
        species = str_replace_all(
            species,
            c(
                "AGLY" = "Aglycines",
                "AGOS" = "Agossypii",
                "APIS" = "Apisum",
                "DNOX" = "Dnoxia",
                "DPLA" = "Dplantaginea",
                "DVIT" = "Dvitifoliae",
                "ELAN" = "Elanigerum",
                "MCER" = "Mcerasi",
                "MPERclone0" = "Mpersicae",
                "PNIG" = "Pnigronervosa",
                "RMAI" = "Rmaidis",
                "RPAD" = "Rpadi",
                "SMIS" = "Smiscanthi"
            )
        )
    ) |>
    mutate(species = factor(species, levels = rev(species_topology))) |>
    ggplot(aes(x = percentage, y = species)) +
    geom_bar(aes(fill = Mode), position = "stack", stat = "identity") +
    facet_wrap(~class, nrow = 1) +
    ggsci::scale_fill_jama() +
    labs(
        x = "Percentage (%)",
        y = ""
    ) +
    theme_bw() +
    theme(axis.text.y = element_blank())

# Combine species tree with percentage plot
final_plot <- wrap_plots(
    p_tree, p_chemosensory_percentage, widths = c(1, 3)
)

final_plot

fig2_A = final_plot





