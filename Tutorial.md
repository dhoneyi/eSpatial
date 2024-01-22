---
title: "R Notebook"
output:
  html_document:
    df_print: paged
---

This is an [R Markdown](http://rmarkdown.rstudio.com) Notebook. When you execute code within the notebook, the results appear beneath the code. 

Try executing this chunk by clicking the *Run* button within the chunk or by placing your cursor inside it and pressing *Ctrl+Shift+Enter*. 

```{r}
obj <- readRDS("../Data/0.obj_melan.rds")


```

```{r}
cols <- c("tumour_1" = "#644498", "tumour_2" = "#70B1D7","Tcell_1"= "#D6ACCF",
         "Tcell_2"="#488CAD","DC/Monocyte"="#7184C1","plasma"="#70B1D7")
DimPlot(obj,reduction = "cod",group.by = "finalType",cols = cols,,pt.size = 1)
```

```{r}
### Calculate spatial specificity scores
clusters <- as.character(unique(obj$finalType))

# Calculate foldchange for RNA
DefaultAssay(obj) <- "RNA"
obj <- NormalizeData(obj)
genes <- rownames(obj)
exp.avg <- data.frame(matrix(nrow = length(genes), ncol = length(clusters)), row.names = genes)
colnames(exp.avg) <- clusters
exp.avg[, clusters] <- sapply(clusters, function(cluster) {
  avg_fc <- FoldChange(obj, ident.1 = cluster)
  return(avg_fc$avg_log2FC)
})
```

```{r}
head(exp.avg)
```


```{r}
# Calculate foldchange for ATAC
DefaultAssay(obj) <- "peaks"
obj <- RunTFIDF(obj)
peaks <- rownames(obj)
peaks.avg <- data.frame(matrix(nrow = length(peaks), ncol = length(clusters)), row.names = peaks)
colnames(peaks.avg) <- clusters
peaks.avg[, clusters] <- t(sapply(clusters, function(cluster) {
  avg_fc <- FoldChange.Seurat(obj, ident.1 = cluster)
  return(avg_fc$avg_log2FC)
}))
```

```{r}
head(peaks.avg)
```

