library(tidyverse)
library(data.table)
library(strataG)
library(nleqslv)
library(ggh4x)
library(rstatix)



# Calcul Kimura distances between TE copies and consensus

# List of enriched transposable elements

list_TE = c("AGLY_TEdenovoGr-B-G3101-Map3", 
                "DPLA_TEdenovoGr-B-G149-Map20_reversed", 
                "DPLA_TEdenovoGr-B-G161-Map14", 
                "DPLA_TEdenovoGr-B-G4075-Map8_reversed",
                "DPLA_TEdenovoGr-B-G4329-Map4",
                "DPLA_TEdenovoGr-B-G77-Map20",
                "DVIT_TEdenovoGr-B-G6314-Map6",
                "DVIT_TEdenovoGr-B-G678-Map3", 
                "DVIT_TEdenovoGr-B-G7829-Map4",
                "DVIT_TEdenovoGr-B-G9099-Map3",
                "MCER_TEdenovoGr-B-G1335-Map3",
                "RPAD_TEdenovoGr-B-G1099-Map3",
                "RPAD_TEdenovoGr-B-G1773-Map4",
                "DVIT_TEdenovoGr-B-G3298-Map9_reversed",
                "DVIT_TEdenovoGr-B-G4429-Map3")

## model "K80" package ape
# K80: The distance derived by Kimura (1980), sometimes referred to as “Kimura's 2-parameters distance”, has the same underlying assumptions than the Jukes–Cantor distance except that two kinds of substitutions are considered: transitions (A <-> G, C <-> T), and transversions (A <-> C, A <-> T, C <-> G, G <-> T). They are assumed to have different probabilities. A transition is the substitution of a purine (C, T) by another one, or the substitution of a pyrimidine (A, G) by another one. A transversion is the substitution of a purine by a pyrimidine, or vice-versa. Both transition and transversion rates are the same for all sites along the DNA sequence. Jin and Nei (1990) modified the Kimura model to allow for variation among sites following a gamma distribution. Like for the Jukes–Cantor model, the gamma parameter must be given by the user. By default, no gamma correction is applied. 

## model "K81" 
# K81: Kimura (1981) generalized his model (Kimura 1980) by assuming different rates for two kinds of transversions: A <-> C and G <-> T on one side, and A <-> T and C <-> G on the other. This is what Kimura called his “three substitution types model” (3ST), and is sometimes referred to as “Kimura's 3-parameters distance”. 

Kimura_res = list()

# A function to calculate Kimura distances between TE associated to chemosensory genes and consensus sequence
Calculate_Kimura_distance = function(TE, substModel = c("K80", "K81"), TE_set = c("chemosensory_associated", "all_TEs")) {
  
  substModel <- match.arg(substModel)
  TE_set <- match.arg(TE_set)
  
  if (TE_set == "chemosensory_associated") {suffix = ".aln"}
  else if (TE_set == "all_TEs") {suffix = "_allTEs.aln"}
  
  sequences_chem = ape::read.dna(paste0("results/",TE, suffix), format = "fasta")
  res = data.frame(TE_consensus = as.character(),
                   TE_copy = as.character(),
                   Ti = as.numeric(),
                   Tv = as.numeric(),
                   Kimura_distance = as.numeric(),
                   TE_set = as.character())
  N = sequences_chem %>% nrow()
  for (j in 1:N) {
    #subset the consensus sequence (the first in the alignment) and one of the others
    sequences_subset = sequences_chem[c(1,j), ]
    
    # Calculate transition/transversion ratio
    TiTv = TiTvRatio(sequences_subset)
    
    # Calculate Kimura distance 
    dist = ape::dist.dna(sequences_subset, model = substModel, pairwise.deletion = TRUE) %>% as.matrix()
    Kimura_dist = dist[1,2]
    
    res_tmp = data.frame(TE_consensus = list_TE[i],
                     TE_copy = labels(sequences_chem)[j],
                     Ti = TiTv[1],
                     Tv = TiTv[2],
                     Kimura_distance = Kimura_dist,
                     TE_set = TE_set)
    res = rbind(res, res_tmp)
  }
  return(res)
}
# A function to calculate transition and transversion rate from distance with Kimura 2 parameters model  
kimura_inverse <- function(K) {
  f <- function(x) {
    p <- x[1] # transition
    q <- x[2] # transversion
    eq1 <- -1/2 * log(1 - 2*p - q) - 1/4 * log(1 - 2*q) - K # eq1 models divergence using the Kimura two-parameter model. it connects p (transition), q (transversion) and K
    eq2 <- p + q - (1/2 * (1 - exp(-4*K)) + 1/2 * (1 - exp(-2*K))) # eq2 ensures that p and q are consistent with the expected values derived from the divergence K
    return(c(eq1, eq2))
  }
  
  res <- nleqslv::nleqslv(c(0.1, 0.1), f)$x # The function uses nleqslv() to solve the nonlinear system. We start with an initial guess (p = 0.1, q = 0.1) and try to adjust p and q to satisfy both equations.
  return(res)
}
# and use these rates to calculate percentage of divergence 
calculate_divergence <- function(K) {
  pq <- kimura_inverse(K)
  p <- pq[1]
  q <- pq[2]
  divergence <- sum(pq)*100
  return(divergence)
}


for (i in 1:length(list_TE)) {
  
  # Start with this TE sequences associated to chemosensory genes
  res1 <- Calculate_Kimura_distance(list_TE[i], substModel = "K80", TE_set = "chemosensory_associated")
  Kimura_res[[length(Kimura_res) + 1]] <- res1
  
  # Then, all sequences of this TE
  res2 <- Calculate_Kimura_distance(list_TE[i], substModel = "K80", TE_set = "all_TEs")
  Kimura_res[[length(Kimura_res) + 1]] <- res2
}

dt_all <- rbindlist(Kimura_res)
dt_all = dt_all[!duplicated(dt_all),]
dt_all = dt_all[! grepl("consensus", TE_copy),]

dt_all <- dt_all %>%
  filter(!is.na(Kimura_distance)) %>% 
  filter(Kimura_distance != "Inf") %>% 
  mutate(divergence = sapply(Kimura_distance, calculate_divergence))

write.table(dt_all, "results/TE_copies_Kimura_distances.csv", row.names = FALSE, sep = "\t")


dt_all = dt_all %>% 
  mutate(species = gsub("_TEdenovoGr-.*", "", TE_consensus)) %>% 
  mutate(TE_set = ifelse(grepl("Or", TE_copy, ignore.case = TRUE), "OR", ifelse(grepl("Gr", TE_copy, ignore.case = TRUE),"GR", TE_set)))


fig_S14 = dt_all %>% 
  ggplot(aes(x = TE_consensus, y = divergence, color = TE_set))+
  geom_boxplot(fill = "grey95", linewidth = 0.5)+
  facet_grid2(~species,  scales = "free")+
  force_panelsizes(cols =  c(1,5,6,1,2))+
  scale_color_manual(values = c("all" = "grey10",
                                "GR" = "dodgerblue3",
                                "OR" = "gold3"))+
  ylab("Kimura distance (%)")+
  xlab("")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1), 
        panel.border = element_rect(colour = "black", fill = NA))

fig3_B = dt_all %>% 
  # filter(Gene_set == "all") %>% 
  ggplot(aes(x = TE_set, y = divergence, color = TE_set))+
  geom_boxplot(fill = "grey95", linewidth = 0.5)+
  scale_color_manual(values = c("all" = "grey10",
                                "GR" = "dodgerblue3",
                                "OR" = "gold3"
  ))+
  ylab("Kimura distance (%)")+
  xlab("")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1), 
        panel.border = element_rect(colour = "black", fill = NA))

ggsave(fig3_B, filename = "results/fig3_B.png", width = 15, height = 30, units = "cm")


# Test statistical differences Kimura distances all TE copies vs GR associated copies vs OR associated copies

# By TE
stat.test <- dt_all %>%
  group_by(TE_consensus) %>%
  do({
    wilcox_test(., divergence ~ TE_set) %>% 
      adjust_pvalue(method = "bonferroni") %>%
      add_significance("p.adj")
  })

# Considering all TEs
stat.test <- dt_all %>%
  do({
    wilcox_test(., divergence ~ TE_set) %>% 
      adjust_pvalue(method = "bonferroni") %>%
      add_significance("p.adj")
  })
