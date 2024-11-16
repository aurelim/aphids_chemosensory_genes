
library("LOLA")
library("GenomicRanges")
library("tidyverse")

## LOLA on aphids

path = "path/to/lola_folder"
sp_list = c("AGOS","AGLY","APIS","DNOX", "DPLA", "DVIT", "ELAN","MCER", "MPER", "PNIG", "RMAI", "RPAD", "SMIS")
gene_set_list = c("_GR", "_OR")
enlargment = c("_2kb", "_10kb")
direction = "enrichment"

for (sp in sp_list) {

  setwd(paste0(path,"/data/data_",sp))
  db_path = paste0(sp, "_TE_collections")
  regionDB = loadRegionDB(db_path, useCache = TRUE)
}

regionDB$regionGRL

lola_res = as.data.frame(matrix(NA, 0,23))
names(lola_res) = c("userSet","dbSet","collection","pValueLog","oddsRatio","support","rnkPV","rnkOR",
                    "rnkSup","maxRnk",  "meanRnk","b","c","d","description","cellType","tissue","antibody",
                    "treatment","dataSource","filename","size", "qvalue")

  
for (sp in sp_list) {
  for (gene_set in gene_set_list) {
    for (enl_size in enlargment) {
      setwd(paste0(path, "data/data_", sp))
      userSets = readBed(paste0(sp, gene_set, enl_size,".sorted.merged.bed"))
      userUniverse = readBed(paste0(sp, enl_size, ".sorted.merged.bed"))
      checkUniverseAppropriateness(userSets, userUniverse)
      
      db_path = paste0(sp, "_TE_collections")
      regionDB = loadRegionDB(db_path, useCache = TRUE)
      
      ## Run analysis
      locResults = runLOLA(userSets, userUniverse, regionDB, cores=2, direction = direction) %>% 
        mutate(Type = str_remove(gene_set, "_"),
               Enlargment_size = str_remove(enl_size, "_"), 
               species = str_remove(collection, "_TE"), 
               filename = str_remove(filename, ".bed"))
      
      locResults$qvalue = p.adjust((10^(-locResults$pValueLog)), method = "bonferroni", n = nrow(locResults))
      
      lola_res = rbind(lola_res, locResults)
      
      writeCombinedEnrichment(locResults, outFolder= paste0(path, "results/",direction,"_july24/merged_", sp, gene_set, enl_size), includeSplits=TRUE)
    }
  }
}

write.table(lola_res,paste0(path, "/results/",direction,"_july24/lola_results.csv"), sep = ";", dec = ".", row.names = FALSE)



lola_res %>% 
  filter(qvalue !=1) %>% 
  ggplot(aes(x = qvalue, fill = userSet))+
  geom_histogram(position = position_dodge())+
  scale_x_log10()+
  facet_wrap(~collection, ncol = 1)+
  geom_vline(xintercept = 0.05)+
  theme_minimal()

lola_res_significatifs = lola_res %>% 
  filter(qvalue <=0.05) %>% 
  select(c(species, Type, Enlargment_size, collection, qvalue, oddsRatio, support, rnkPV, rnkOR, rnkSup, maxRnk, meanRnk, b, c, d, filename, size))

write.table(lola_res_significatifs, file = paste0(path, "/results/",direction,"_july24/lola_signif.csv"), sep = ";", dec = ".", row.names = FALSE)
            

