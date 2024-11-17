library(tidyverse)
library(data.table)
library(seqinr)
library(reshape2)
library(plotly)
library(ggh4x)
library(RColorBrewer)

path = "/home/aurelie/Dropbox/Post-Doc_ECLECTIC/8_Aphids_comparative_genomic/tegrip_aphids/"

list_TE_lola = c("AGLY_TEdenovoGr-B-G3101-Map3",
                 "DPLA_TEdenovoGr-B-G149-Map20_reversed",
                 "DPLA_TEdenovoGr-B-G161-Map14",
                 "DPLA_TEdenovoGr-B-G4075-Map8_reversed",
                 "DPLA_TEdenovoGr-B-G4329-Map4",
                 "DPLA_TEdenovoGr-B-G77-Map20",
                 "DVIT_TEdenovoGr-B-G6314-Map6",
                 "DVIT_TEdenovoGr-B-G3298-Map9_reversed",
                 "DVIT_TEdenovoGr-B-G4429-Map3",
                 "DVIT_TEdenovoGr-B-G678-Map3",
                 "DVIT_TEdenovoGr-B-G7829-Map4",
                 "DVIT_TEdenovoGr-B-G9099-Map3",
                 "MCER_TEdenovoGr-B-G1335-Map3",
                 "RPAD_TEdenovoGr-B-G1099-Map3",
                 "RPAD_TEdenovoGr-B-G1773-Map4")

sp_list = c("AGLY", "AGOS", "DNOX", "DPLA", "DVIT", "ELAN", "MCER", "MPER", "PNIG", "RMAI", "RPAD", "SMIS")

colors_list = brewer.pal(n = 6, name = "Spectral")


# A function inspired by TEgrip algorithm to detect all relationships between genes of interest and enriched TE detected using LOLA
calcul_distance = function(TE, Gene_family) {
  res <- data.frame(Chr = character(), TE_id = character(), TE_strand = character(), TE_start = numeric(), TE_end = numeric(),
                    Gene_id = character(), Gene_strand = character(), Gene_start = numeric(), Gene_end = numeric(),
                    distance1 = numeric(), distance2 = numeric(), distance3 = numeric(), distance4 = numeric(),
                    relationship = character(), Gene_fam = character(), stringsAsFactors = FALSE)
  
  species = gsub("_TEdenovoGr.*", "", TE)
  TE_input = fread(paste0(path, "data/",species,"_TE_annotation.tsv")) %>% 
    select(c(seqid, start, end, strand, class, TEname)) %>% filter(TEname == TE)
  
  gene_input = fread(paste0(path, "data/",species,"_",Gene_family,"gene_annotation.tsv")) %>% select(c(chromosome, start, end, strand, ID))
  liste_genes = unique(gene_input$ID)
  
  for (gene in liste_genes) {
    
    gene_info = filter(gene_input, ID == gene)
    
    chromosome = gene_info$chromosome
    start = gene_info$start
    end = gene_info$end
    strand = gene_info$strand
    
    TE_in_same_chr = filter(TE_input, seqid == chromosome)
    
    if(nrow(TE_in_same_chr) > 0) {
      tmp <- data.frame(
        Chr = chromosome,
        TE_id = TE,
        TE_strand = TE_in_same_chr$strand,
        TE_start = TE_in_same_chr$start,
        TE_end = TE_in_same_chr$end,
        Gene_start = start,
        Gene_end = end,
        Gene_strand = strand,
        distance1 = TE_in_same_chr$start - end,
        distance2 = start - TE_in_same_chr$end,
        distance3 = end - TE_in_same_chr$end,
        distance4 = start - TE_in_same_chr$start,
        Gene_id = gene,
        relationship = "",
        Gene_fam = Gene_family
      )
      
      tmp <- tmp %>%
        rowwise() %>%
        mutate(relationship = case_when(
          distance1 < 0 & distance2 < 0 & distance3 >= 0 & distance4 <= 0 ~ "subset",
          distance1 < 0 & distance2 < 0 & distance3 <= 0 & distance4 >= 0 ~ "superset",
          distance1 < 0 & distance2 < 0 & distance3 < 0 & distance4 < 0 ~ "downstream_overlap",
          distance1 < 0 & distance2 < 0 & distance3 > 0 & distance4 > 0 ~ "upstream_overlap",
          distance1 > 0 & distance2 < 0 & distance3 < 0 & distance4 < 0 ~ "downstream",
          distance1 < 0 & distance2 > 0 & distance3 > 0 & distance4 > 0 ~ "upstream",
          TRUE ~ ""
        )) %>%
        ungroup() %>%
        mutate(relationship = ifelse((relationship == "upstream" | relationship == "downstream") & 
                                       pmin(abs(distance1), abs(distance2), abs(distance3), abs(distance4)) > 10000, "", relationship))
      
      tmp = tmp %>% filter(relationship != "")
      
      # Trouver la ligne avec la plus petite valeur parmi les colonnes sélectionnées
      min_row_index <- which.min(apply(tmp[,c("distance1", "distance2", "distance3", "distance4")], 1, 
                                       function(x) min(abs(x))))
      
      # Extraire la ligne correspondante
      result <- tmp[min_row_index, ]
      
      res <- bind_rows(res, tmp)
    } else {res = res}
  }
  return(res)
}


results_list <- list()

for (TE in list_TE_lola) {

  res_OR = calcul_distance(TE, "OR") %>% 
    filter(relationship != "") %>% 
    select(-c(distance1, distance2, distance3, distance4))
  
  results_list[[length(results_list) + 1]] <- res_OR
  
  res_GR = calcul_distance(TE, "GR") %>% 
    filter(relationship != "") %>% 
    select(-c(distance1, distance2, distance3, distance4))
  
  results_list[[length(results_list) + 1]] <- res_GR
}

dt_all <- rbindlist(results_list)
dim(dt_all)
head(dt_all)
write.table(dt_all, file = paste0(path, "results/Enriched_te.csv"), row.names = FALSE, sep = "\t")


# Count gene number by TE and TE number by gene 

Count_genes = dt_all %>% filter(TE_id %in% list_TE_lola) %>% 
  mutate(species = str_remove_all(TE_id, "_TEdenovoGr.*")) %>% 
  group_by(TE_id, Gene_fam, relationship, species) %>% 
  dplyr::count() %>% 
  select(c(species, TE_id, Gene_fam, relationship, n))
write.table(Count_genes, file = paste0(path, "results/Count_genes.csv"), row.names = FALSE, sep = "\t")


Count_TE = dt_all %>% filter(TE_id %in% list_TE_lola) %>% 
  mutate(species = str_remove_all(TE_id, "_TEdenovoGr.*")) %>% 
  group_by(Gene_id, Gene_fam, relationship, species) %>% 
  dplyr::count()
write.table(Count_TE, file = paste0(path, "results/Count_TE.csv"), row.names = FALSE, sep = "\t")
  


fig3_A = Count_genes %>% 
  ggplot(aes(x = n, y = TE_id, fill = relationship))+
  geom_bar(stat = "identity")+
  facet_grid(species~Gene_fam, scales = "free_y")+
  force_panelsizes(rows = c(1,5,6,1,2))+
  scale_fill_manual(values = c("downstream" = colors_list[1],
                               "downstream_overlap" = colors_list[2],
                               "subset" = colors_list[3],
                               "superset" = colors_list[4],
                               "upstream" = colors_list[5],
                               "upstream_overlap" = colors_list[6]))+
  scale_y_discrete(labels = 
                     c("AGLY_TEdenovoGr-B-G3101-Map3" = "AGLY-G3101-Map3 | TIR hAT",
                       "DPLA_TEdenovoGr-B-G149-Map20_reversed" = "DPLA-G149-Map20 | TIR Tc1-Mariner",
                       "DPLA_TEdenovoGr-B-G161-Map14" = "DPLA-G161-Map14 | class II",
                       "DPLA_TEdenovoGr-B-G4075-Map8_reversed" = "DPLA-G4075-Map8 | SINE",
                       "DPLA_TEdenovoGr-B-G4329-Map4" = "DPLA-G4329-Map4 | TIR hAT" ,
                       "DPLA_TEdenovoGr-B-G77-Map20" = "DPLA-G77-Map20 | class II",
                       "DVIT_TEdenovoGr-B-G6314-Map6" = "DVIT-G6314-Map6 | class II",
                       "DVIT_TEdenovoGr-B-G3298-Map9_reversed" = "DVIT-G3298-Map9 | TIR hAT",
                       "DVIT_TEdenovoGr-B-G4429-Map3" = "DVIT-G4429-Map3 | TIR hAT",
                       "DVIT_TEdenovoGr-B-G678-Map3" = "DVIT-G678-Map3 | TIR Mutator",
                       "DVIT_TEdenovoGr-B-G7829-Map4" = "DVIT-G7829-Map4 | TIR hAT",
                       "DVIT_TEdenovoGr-B-G9099-Map3" = "DVIT-G9099-Map3 | TIR hAT",
                       "MCER_TEdenovoGr-B-G1335-Map3" = "MCER-G1335-Map3 | SINE",
                       "RPAD_TEdenovoGr-B-G1099-Map3" = "RPAD-G1099-Map3 | TIR Tc1-Mariner",
                       "RPAD_TEdenovoGr-B-G1773-Map4" = "RPAD-G1773-Map4 | unclassified")) +
  
  ylab("Enriched TEs")+
  xlab("gene count")+
  theme_linedraw()+
  theme(strip.background = element_rect(fill = "grey80", colour = "grey30"),
        strip.text = element_text(color = "black"))

ggsave(fig3_A, filename = "results/fig3_A.png", width = 15, height = 30, units = "cm")


# Percentage of genes with associated TEs
TEenriched_genes = unique(dt_all$Gene_id)

Chemgenes_annotation = do.call(rbind, lapply(seq_along(sp_list), function(x){
  sp = sp_list[x]

  ORgenes = fread(paste0(path, "data/",sp,"_ORgene_annotation.tsv")) %>% 
    mutate(species = sp) %>% 
    mutate(gene_class = "OR")
  GRgenes = fread(paste0(path, "data/",sp,"_GRgene_annotation.tsv")) %>% 
    mutate(species = sp) %>% 
    mutate(gene_class = "GR")
  tmp = rbind(ORgenes, GRgenes)
  return(tmp)
})) %>% 
  mutate(TE_enriched_gene = ifelse(ID %in% TEenriched_genes, "yes", "no"))

Chemgenes_annotation %>% 
  group_by(species, TE_enriched_gene, gene_class) %>% 
  tally() %>% 
  print(n = 32)
