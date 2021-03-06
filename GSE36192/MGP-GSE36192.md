---
title: "MGP GSE36192"
author: "Keon Arbabi"
date: "06/05/2021"
output: 
  html_document: 
    keep_md: yes
---

# Knitr settings 



# Load packages 


```r
# # Installs packages from Bioconductor
# BiocManager::install("GEOquery")
# BiocManager::install("edgeR")
# BiocManager::install("annotationTools")
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
```

# Load data   

Data processing: "Illumina Genome Studio Gene Expression Module v3.2.7, cubic spline normalization"


```r
# load series and platform data from GEO
gset = getGEO("GSE36192", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx = grep("GPL6947", attr(gset, "names")) else idx = 1
gset = gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) = make.names(fvarLabels(gset))

# group membership for all samples
gsms = paste0("11111111111111111111111111111111111111111111111111",
        "11111111111111111111111111111111111111111111111111",
        "11111111111111111111111111111111111111111111111111",
        "11111111111111111111111111111111111111111111111111",
        "11111111111111111111111111111111111111111111111111",
        "11111111111111111111111111111111111111111111111111",
        "11111111111111111111111111111111111111111111111111",
        "11111111111111111111111111111111111111111111111111",
        "11111111111111111111111111111111111111111111111111",
        "11111100000000000000000000000000000000000000000000",
        "00000000000000000000000000000000000000000000000000",
        "00000000000000000000000000000000000000000000000000",
        "00000000000000000000000000000000000000000000000000",
        "00000000000000000000000000000000000000000000000000",
        "00000000000000000000000000000000000000000000000000",
        "00000000000000000000000000000000000000000000000000",
        "00000000000000000000000000000000000000000000000000",
        "00000000000000000000000000000000000000000000000000",
        "00000000000")
sml = strsplit(gsms, split="")[[1]]

pData(gset)$data_processing[1]
```

```
## [1] "Illumina Genome Studio Gene Expression Module v3.2.7, cubic spline normalization"
```

```r
# log2 transformation
ex = exprs(gset)
qx = as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC = (qx[5] > 100) ||
          (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { ex[which(ex <= 0)] = NaN
  exprs(gset) = log2(ex) }

# assign samples to groups and set up design matrix
gs = factor(sml)
groups = make.names(c("frontal cortex","cerebellum"))
levels(gs) = groups
gset$group = gs
design = model.matrix(~group + 0, gset)
colnames(design) = levels(gs)

fit = lmFit(gset, design)  # fit linear model

# set up contrasts of interest and recalculate model coefficients
cts = paste(groups[1], groups[2], sep="-")
cont.matrix = makeContrasts(contrasts=cts, levels=design)
fit2 = contrasts.fit(fit, cont.matrix)

# compute statistics and table of top significant genes
fit2 = eBayes(fit2, 0.01)
tT = topTable(fit2, adjust="fdr", sort.by="B", number=250)

tT = subset(tT, select=c("ID","adj.P.Val","P.Value","t","B","logFC","Gene.symbol","Gene.title"))
#write.table(tT, file=stdout(), row.names=F, sep="\t")

# Visualize and quality control test results.
# Build histogram of P-values for all genes. Normal test
# assumption is that most genes are not differentially expressed.
tT2 = topTable(fit2, adjust="fdr", sort.by="B", number=Inf)
hist(tT2$adj.P.Val, col = "grey", border = "white", xlab = "P-adj",
  ylab = "Number of genes", main = "P-adj value distribution")
```

![](MGP-GSE36192_files/figure-html/unnamed-chunk-2-1.png)<!-- -->

```r
# summarize test results as "up", "down" or "not expressed"
dT = decideTests(fit2, adjust.method="fdr", p.value=0.05)

# Venn diagram of results
vennDiagram(dT, circle.col=palette())
```

![](MGP-GSE36192_files/figure-html/unnamed-chunk-2-2.png)<!-- -->

```r
# create Q-Q plot for t-statistic
t.good = which(!is.na(fit2$F)) # filter out bad probes
qqt(fit2$t[t.good], fit2$df.total[t.good], main="Moderated t statistic")
```

![](MGP-GSE36192_files/figure-html/unnamed-chunk-2-3.png)<!-- -->

```r
# volcano plot (log P-value vs log fold change)
colnames(fit2) # list contrast names
```

```
## [1] "frontal.cortex-cerebellum"
```

```r
ct = 1        # choose contrast of interest
volcanoplot(fit2, coef=ct, main=colnames(fit2)[ct], pch=20,
  highlight=length(which(dT[,ct]!=0)), names=rep('+', nrow(fit2)))
```

![](MGP-GSE36192_files/figure-html/unnamed-chunk-2-4.png)<!-- -->

```r
# MD plot (log fold change vs mean log expression)
# highlight statistically significant (p-adj < 0.05) probes
plotMD(fit2, column=ct, status=dT[,ct], legend=F, pch=20, cex=1)
abline(h=0)
```

![](MGP-GSE36192_files/figure-html/unnamed-chunk-2-5.png)<!-- -->

```r
################################################################
# General expression data analysis
ex = exprs(gset)

# box-and-whisker plot
dev.new(width=3+ncol(gset)/6, height=5)
ord = order(gs)  # order samples by group
ord = ord[c(1:100,(length(ord)-100):(length(ord)))]
palette(c("#1B9E77", "#7570B3", "#E7298A", "#E6AB02", "#D95F02",
          "#66A61E", "#A6761D", "#B32424", "#B324B3", "#666666"))
par(mar=c(7,4,2,1))
title = paste ("GSE36192", "/", annotation(gset), sep ="")
boxplot(ex[,ord], boxwex=0.6, notch=T, main=title, outline=FALSE, las=2, col=gs[ord])
legend("topleft", groups, fill=palette(), bty="n")
dev.off()
```

```
## png 
##   2
```

```r
# expression value distribution
par(mar=c(4,4,2,1))
title = paste ("GSE36192", "/", annotation(gset), " value distribution", sep ="")
plotDensities(ex, group=gs, main=title, legend ="topright")
```

![](MGP-GSE36192_files/figure-html/unnamed-chunk-2-6.png)<!-- -->

```r
# UMAP plot (dimensionality reduction)
ex = na.omit(ex) # eliminate rows with NAs
ex = ex[!duplicated(ex), ]  # remove duplicates
ump = umap(t(ex), n_neighbors = 15, random_state = 123)
par(mar=c(3,3,2,6), xpd=TRUE)
plot(ump$layout, main="UMAP plot, nbrs=15", xlab="", ylab="", col=gs, pch=20, cex=1.5)
legend("topright", inset=c(-0.15,0), legend=levels(gs), pch=20, col=1:nlevels(gs), title="Group", pt.cex=1.5)
```

![](MGP-GSE36192_files/figure-html/unnamed-chunk-2-7.png)<!-- -->

```r
# mean-variance trend, helps to see if precision weights are needed
plotSA(fit2, main="Mean variance trend, GSE36192")
```

![](MGP-GSE36192_files/figure-html/unnamed-chunk-2-8.png)<!-- -->

# Build data frames


```r
# Build sample metadata
data_meta = pData(gset)
data_meta = data.frame(geo_accession = gset$geo_accession,
                       age_years = gset$`age (y):ch1`,
                       sex = gset$`gender:ch1`,
                       region = gset$`tissue:ch1`,
                       pmi = gset$`pmi (hr):ch1`
                       ) 
data_meta %<>% mutate_at(c("age_years","pmi"), ~as.numeric(as.character(.)))

# View distribution of ages 
ggplot(data_meta %>% filter(region == "frontal cortex"), aes(age_years)) + 
  geom_histogram(bins = 50) +
  theme_bw()
```

![](MGP-GSE36192_files/figure-html/unnamed-chunk-3-1.png)<!-- -->

```r
paste("frontal cortex individuals:", nrow(data_meta %>% filter(region == "frontal cortex")))
```

```
## [1] "frontal cortex individuals: 455"
```

```r
paste("sample age mean: ", mean(unlist(data_meta %>% pull(age_years))))
```

```
## [1] "sample age mean:  48.8173874862788"
```

```r
paste("sample age sd: ", sd(unlist(data_meta %>% pull(age_years))))
```

```
## [1] "sample age sd:  25.6334661357824"
```

```r
# Grab probe annotations
if(is(gset,"ExpressionSet")) {
	y = get("exprs",env=gset@assayData)
	if(length(gset@featureData@data)) probe_ann = gset@featureData@data
	if(!is.null(rownames(y))) {
		if(is.null(probe_ann))
			probe_ann = data.frame(ID=I(rownames(y)))
		else
			probe_ann$ID = rownames(y)
		}}
probe_ann %<>% select(ID, Gene.symbol)

# Combined with probe list to get gene symbols
data_exp = ex %>% as.data.frame() %>% rownames_to_column(var = "ID")
data_exp = inner_join(probe_ann, data_exp, by = "ID")
# Remove blank probes 
data_exp = data_exp[-which(data_exp$Gene.symbol == ""),]

# Collapse multiple microarray probes for a single gene and then merge the data by gene identifier; choosing the probe with the highest 
# average expression leads to best between-study consistency 
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3166942/
rownames(data_exp) = data_exp$ID
select.rows = collapseRows(datET = data_exp[,-c(1:2)],
                    rowGroup = data_exp$Gene.symbol,
                    rowID = data_exp$ID,
                    method = "maxRowVariance")
rownames(data_exp) = NULL
data_exp %<>% filter(ID %in% select.rows$group2row[,2]) %>%
  dplyr::select(-ID) 

# # Calculate median expression
# medExp = data_exp%>%
#     sepExpr() %>% {.[[2]]} %>%
#     unlist %>% median
# # Remove genes with a median expression below the median expression of the dataset
# data_exp = mostVariable(data_exp, genes = "gene_assignment", threshold = medExp, threshFun = median)
# data_exp = data_exp[complete.cases(data_exp), ]
# rownames(data_exp) = NULL

# Merge gene expression and meta dataframes
data_comb = data_exp %>%   
  column_to_rownames(var = "Gene.symbol") %>%
  t() %>% as.data.frame() %>%
  rownames_to_column(var = "geo_accession")

data_comb = inner_join(data_meta, data_comb, by = 'geo_accession')
write_rds(data_comb, file = here("output", "GSE36192_out.rds"))
```

# Cell type proportion estimation


```r
# Filter data
data_working = data_comb %>% filter(region == "frontal cortex") %>%
  column_to_rownames(var = "geo_accession") %>%
  select(-colnames(data_meta[,-1])) %>%
  t() %>% as.data.frame() %>%
  rownames_to_column(var = "Gene.Symbol")

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

paste("marker matches in data: ", length(intersect(unlist(data_working$Gene.Symbol), unlist(marker_data$gene))), "/",
      nrow(marker_data))
```

```
## [1] "marker matches in data:  1088 / 1357"
```

```r
# Get vector of unique cell types 
(cell_types = marker_data$subclass %>% unique())
```

```
##  [1] "VLMC"             "Oligodendrocyte"  "Endothelial"      "Inh_PAX6"        
##  [5] "Exc_L5/6_IT_Car3" "Microglia"        "Pericyte"         "Astrocyte"       
##  [9] "Exc_IT"           "Exc_L6_CT"        "OPC"              "Exc_L5_ET"       
## [13] "Inh_PVALB"        "Inh_VIP"          "Exc_L5/6_NP"      "Inh_LAMP5"       
## [17] "Inh_SST"          "Exc_L6b"          "Exc_L4_IT"
```

```r
# Organize markers into a list 
marker_list = lapply(cell_types, function(cell_type){
  return(marker_data %>% filter(subclass == cell_type) %>% pull(gene) %>% unlist())
  })
names(marker_list) = cell_types

# Run MGP analysis
estimations =  mgpEstimate(
  exprData = data_working,
  genes = marker_list,
  geneColName = 'Gene.Symbol',
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
  rownames_to_column(var = "geo_accession")
mgp_df = inner_join(data_meta, mgp_estimates, by = "geo_accession") %>%
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
plot_grid(plotlist = plot_list, nrow = 1)
```

![](MGP-GSE36192_files/figure-html/unnamed-chunk-4-1.png)<!-- -->

# Linear models


```r
# Linear models where cell_type_prop ~ sex + pmi + age_years`
lm_df = mgp_df %>%
  group_by(cell_type) %>%
  do(tidy(lm(scale(cell_proportion) ~ sex + scale(pmi) + scale(age_years),  data = .))) %>%
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

![](MGP-GSE36192_files/figure-html/unnamed-chunk-5-1.png)<!-- -->


