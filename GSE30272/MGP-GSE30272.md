---
title: "MGP GSE30272"
author: "Keon Arbabi"
date: "16/04/2021"
output: 
  html_document: 
    keep_md: yes
---

# Knitr settings 



# Load packages 


```r
## Installs packages from Bioconductor 
# BiocManager::install("GEOquery")
# BiocManager::install("edgeR")
# install.packages("cli") 
# install.packages("processx")
# install.packages("rlang")
# devtools::install_github('oganm/markerGeneProfile', force = T)

# Load packages 
library(GEOquery) 
library(limma)
library(umap)
library(tidyverse) 
library(magrittr) 
library(ggpubr) 
library(splitstackshape) 
library(WGCNA) 
library(here)
library(edgeR)
library(markerGeneProfile)
library(matrixStats)
library(cowplot)
library(broom)
library(annotationTools)
library(maptools)
library(RCurl)
```

# Load and process metadata


```r
# Download from GSE
gset = getGEO(GEO = "GSE30272", filename = NULL, destdir = "./", GSElimits = NULL, GSEMatrix = TRUE, AnnotGPL = FALSE, getGPL = FALSE)
# Unlist to get ExpressionSet object
gset = gset[[1]]
pData(gset)$data_processing[1]
```

```
## [1] "Two color intensity data was imported, background subtracted, low intensities dropped, intensities converted to log2(sample/ref), ratios loess normalized, outliers dropped (>6 MADs), missing data from low intesity cut-off were imputed. Only the 30,176 probes (out of the total 49,152 on GPL4611) that passed all quality control and other filters were used. No \"cleaning\" procedure has been applied to these data - see \"Overall design\" and \"Supplemental Files\" for additional data."
```

```r
# Get identifiers from pData
data_meta = data.frame(geo_accession = gset$geo_accession,
                       sample_id = gset$title,
                       age_years = gset$characteristics_ch1.1,
                       sex = gset$characteristics_ch1.2,
                       pmi = gset$characteristics_ch1.4,
                       rin = gset$characteristics_ch1.6,
                       ph = gset$`ph:ch1`,
                       batch =  gset$characteristics_ch1)
# Clean up variables a little 
data_meta = cSplit(indt = data_meta, splitCols = c("age_years","sex","pmi","rin","batch"), 
                  sep = ":", direction = "wide", stripWhite = TRUE)
data_meta = data_meta[,-c(4,6,8,10,12)]
colnames(data_meta) = c("geo_accession","sample_id","ph","age_years","sex","pmi","rin","batch")
data_meta %<>% relocate(ph, .after = rin)

# Get sample ids for prenatal samples
postnatal_samples = data_meta %>% filter(age_years > 15) %>% pull(sample_id)

# View distribution of ages 
ggplot(data_meta, aes(age_years)) + 
  geom_histogram(bins = 100) +
  theme_bw()
```

![](MGP-GSE30272_files/figure-html/unnamed-chunk-2-1.png)<!-- -->

```r
ggplot(data_meta %>% filter(sample_id %in% postnatal_samples), aes(age_years)) + 
  geom_histogram(bins = 80) +
  theme_bw()
```

![](MGP-GSE30272_files/figure-html/unnamed-chunk-2-2.png)<!-- -->

```r
paste("sample age mean: ", mean(unlist(data_meta %>% filter(age_years > 15) %>% pull(age_years))))
```

```
## [1] "sample age mean:  39.4626948121858"
```

```r
paste("sample age sd: ", sd(unlist(data_meta %>% filter(age_years > 15) %>% pull(age_years))))
```

```
## [1] "sample age sd:  16.8421875767115"
```

# Load and process expression data

This is a version of the data tuned by the original authors to identifying canonical patterns of gene expression across the lifespan (at the expense of individual variation). It expands on previous modeling (Colantuoni 2011, PMID: 22031444) of age patterns across the lifespan in several key ways (manuscript in preparation): 1) We applied splines to capture non-linear gene expression effects while ensuring patterns of gene expression are continuous across the lifespan. The previous analysis used age by decade interaction terms, which are not necessarily continuous. 2) We estimated and adjusted for a much higher number of SVs. The previous analysis used only 2 SVs, here we allowed SVA to automatically determine this number: 31 SVs were used. This much increased "cleaning" further tuned this dataset to age effects. Hence, this newly processed data should only be used for the estimation of canonical, mean patterns of expression across the lifetime. 3) We regressed out SVs while allowing the effects of age and mean gene expression (the intercept) to remain in the data. Previously, SVs were regressed out while ignoring possible correlation between SVs and age, potentially obscuring some age effects. Specifically, using SVA, we employ a 2nd degree basis spline with knots at birth, 1, 10, 20, and 50 years [8 degrees of freedom], i.e. a curve fit to expression across age within each age range between these knots. Each model also allowed an offset at birth, because there were no samples in the third trimester of fetal life


```r
# Cleaned data provided by the authors https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE30272
data_exp = read.csv(file = here("input", "GSE30272_ExprsMtxCleanedN269_31SVN.csv"))
colnames(data_exp)[1] = "probe_id"
# Metadata annotating probe ids to genes https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL4611
annot = read.csv(file = here("input", "GSE30272_annot.csv"))

data_exp %<>% as.data.frame() %>%
  mutate(probe_id = as.character(probe_id)) %>%
  mutate(gene_name = annot$Gene_Symbol[match(probe_id, annot$ID)]) %>%
  dplyr::select(gene_name, probe_id, everything())

# Remove rows with missing data 
data_exp = data_exp[-which(data_exp == ""),]
data_exp = data_exp[-which(data_exp == "##noname##"),]
# Double-check
data_exp = data_exp[complete.cases(data_exp),]
data_exp %<>% mutate_if(is.factor, ~as.numeric(as.character(.)))

# Collapse multiple microarray probes for a single gene and then merge the data by gene identifier; choosing the probe with the highest 
# average expression leads to best between-study consistency 
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3166942/
rownames(data_exp) = data_exp$probe_id
select.rows = collapseRows(datET = data_exp[,-c(1:2)],
                    rowGroup = data_exp$gene_name,
                    rowID = data_exp$probe_id,
                    method = "MaxMean")
rownames(data_exp) = NULL
data_exp %<>% filter(probe_id %in% select.rows$group2row[,2]) %>%
  dplyr::select(-probe_id) 

# box-and-whisker plot
title = paste ("GSE30272", "/", annotation(gset), sep ="")
boxplot(data_exp[-1], boxwex=0.7, notch=T, main=title, outline=FALSE, las=2)
```

![](MGP-GSE30272_files/figure-html/unnamed-chunk-3-1.png)<!-- -->

```r
# expression value distribution plot
par(mar=c(4,4,2,1))
title = paste ("GSE30272", "/", annotation(gset), " value distribution", sep ="")
plotDensities(data_exp[-1], main=title, legend=F)
```

![](MGP-GSE30272_files/figure-html/unnamed-chunk-3-2.png)<!-- -->

```r
# mean-variance trend
plotSA(lmFit(data_exp[-1]), main="Mean variance trend, GSE30272")
```

![](MGP-GSE30272_files/figure-html/unnamed-chunk-3-3.png)<!-- -->

```r
# UMAP plot (multi-dimensional scaling)
ump = umap(t(data_exp[-1]), n_neighbors = 15, random_state = 123)
plot(ump$layout, main="UMAP plot, nbrs=15", xlab="", ylab="", pch=20, cex=1.5)
```

![](MGP-GSE30272_files/figure-html/unnamed-chunk-3-4.png)<!-- -->

```r
# Merge gene expression and meta dataframes
data_comb = data_exp %>%   
  column_to_rownames(var = "gene_name") %>%
  t() %>% as.data.frame() %>%
  rownames_to_column(var = "sample_id")

data_comb = inner_join(data_meta, data_comb, by = 'sample_id')
write_rds(data_comb, file = here("output", "GSE30272_out.rds"))
```

# Cell type proportion estimation


```r
# Load in marker genes
marker_data = read.csv("https://github.com/sonnyc247/MarkerSelection/raw/master/Data/Outputs/CSVs_and_Tables/Markers/MTG_and_CgG_lfct2/new_MTGnCgG_lfct2.5_results.csv")

# Manually add inhibitory/excitatory suffix to subclass labels 
marker_data %<>%
  mutate(class = case_when(
    subclass == "LAMP5" ~ "Inh",
    subclass == "L5 ET" ~ "Exc",
    subclass == "PAX6" ~ "Inh",
    subclass == "L5/6 IT Car3" ~ "Exc",
    subclass == "VIP" ~ "Inh",
    subclass == "L6 CT" ~ "Exc",
    subclass == "IT"  ~ "Exc",
    subclass == "L6b" ~ "Exc",
    subclass == "PVALB" ~ "Inh",
    subclass == "L5/6 NP" ~ "Exc",
    subclass == "SST" ~ "Inh",
    subclass == "L4 IT" ~ "Exc"
  )) %>% 
  relocate(class, .before = "subclass") %>%
  unite(subclass, c(class, subclass), sep = "_", remove = F, na.rm = T)
marker_data$subclass = gsub(" ", "_", marker_data$subclass)
marker_data$class[is.na(marker_data$class)] = "NonN"

paste("marker matches in data: ", length(intersect(unlist(data_exp$gene_name), unlist(marker_data$gene))), "/",
      nrow(marker_data))
```

```
## [1] "marker matches in data:  944 / 1357"
```

```r
# Get vector of unique cell types 
cell_types = marker_data$subclass %>% unique()
# Organize markers into a list 
marker_list = lapply(cell_types, function(cell_type){
  return(marker_data %>% filter(subclass == cell_type) %>% pull(gene) %>% unlist())
  })
names(marker_list) = cell_types

# Run MGP analysis
estimations =  mgpEstimate(
  exprData = data_exp[c("gene_name", postnatal_samples)],
  genes = marker_list,
  geneColName = 'gene_name',
  outlierSampleRemove = FALSE, # should outlier samples removed. This is done using boxplot stats
  geneTransform = NULL, # this is the default option for geneTransform
  groups = NULL, # if there are experimental groups provide them here. if not desired set to NULL
  seekConsensus = FALSE, # ensures gene rotations are positive in both of the groups
  removeMinority = TRUE)
# Lost cell types 
setdiff(cell_types, names(estimations$estimates))
```

```
## character(0)
```

```r
# Merge cell type proportions with sample metadata
mgp_estimates = as.data.frame(estimations$estimates) %>%
  rownames_to_column(var = "sample_id")
mgp_df = inner_join(data_meta, mgp_estimates, by = "sample_id") %>%
  pivot_longer(-colnames(data_meta),
               names_to = "cell_type",
               values_to = "cell_proportion")

plot_genes = c("Exc_IT","Inh_SST","Oligodendrocyte")
plot_list = list()

for(i in 1:length(plot_genes)){
  plot_list[[i]] = ggplot(
    mgp_df %>% filter(cell_type == plot_genes[i]),
    aes(x = age_years, y = cell_proportion)) +
    geom_smooth(method = "lm", se = F) + 
    geom_point() +
    ylab(paste(plot_genes[i], " MGP")) + xlab ("sample age (years)") +
    theme_bw()
}
plot_grid(plotlist = plot_list, nrow =1)
```

![](MGP-GSE30272_files/figure-html/unnamed-chunk-4-1.png)<!-- -->

# linear models


```r
# Linear models where cell_type_prop ~ sex + rin + pmi + age_years`
lm_df = mgp_df %>%
  group_by(cell_type) %>%
  do(tidy(lm(scale(cell_proportion) ~ sex + scale(rin) + scale(pmi) + scale(age_years),  data = .))) %>%
  ungroup() %>%
  mutate(padj = p.adjust(`p.value`, method = 'BH')) %>%
  mutate(term = recode(term, 
                       `(Intercept)` = "Intercept", 
                       `sexM` = "sex:M",
                       `sexF` = "sex:F",
                       `scale(rin)` = "rin",
                       `scale(pmi)` = "pmi",
                       `scale(age_years)` = "age_years")) %>%
  mutate(class = case_when(
    str_detect(cell_type, "Inh") ~ "Inhibitory",
    str_detect(cell_type, "Exc") ~ "Excitatory",
    TRUE ~ "Non-Neuronal"
  ))

# Beta coeffs per cell type for phenotype and age effects
beta_plot = lm_df %>% 
  filter(term %in% 'age_years') %>% 
  mutate(cell_type = fct_reorder(cell_type, estimate)) %>% 
  ggplot(aes(x = cell_type, y = estimate)) + 
  geom_hline(yintercept = 0) + 
  geom_bar(stat = "identity") + 
  geom_errorbar(aes(ymin = estimate - std.error, ymax = estimate + std.error)) + 
  ylab('Std. Beta coeff.') + 
  xlab('Cell type proportions') + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  facet_wrap(~class, drop = T, scale = "free")

beta_plot
```

![](MGP-GSE30272_files/figure-html/unnamed-chunk-5-1.png)<!-- -->

