rm(list = ls())
# setwd("D:\\300_学习\\340_博士后\\Deer\\deer")
setwd("/public6/yuanming/paper/Deer/deer")

suppressMessages(suppressWarnings(library(devtools)))
suppressMessages(suppressWarnings(library(dplyr)))
load_all()




### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### mmu
### 

wholeFilesPath <- c(
    "/public6/yuanming/paper/publicData/BM1_liver_0d_rsem/BM1_liver_0d_rsem.genes.results",
    "/public6/yuanming/paper/publicData/BM2_liver_0d_rsem/BM2_liver_0d_rsem.genes.results",
    "/public6/yuanming/paper/publicData/BM3_liver_0d_rsem/BM3_liver_0d_rsem.genes.results",
    "/public6/yuanming/paper/publicData/BM1_liver_7d_rsem/BM1_liver_7d_rsem.genes.results",
    "/public6/yuanming/paper/publicData/BM2_liver_7d_rsem/BM2_liver_7d_rsem.genes.results",
    "/public6/yuanming/paper/publicData/BM3_liver_7d_rsem/BM3_liver_7d_rsem.genes.results",
    "/public6/yuanming/paper/publicData/BM1_liver_10d_rsem/BM1_liver_10d_rsem.genes.results",
    "/public6/yuanming/paper/publicData/BM2_liver_10d_rsem/BM2_liver_10d_rsem.genes.results",
    "/public6/yuanming/paper/publicData/BM3_liver_10d_rsem/BM3_liver_10d_rsem.genes.results",
    "/public6/yuanming/paper/publicData/BM1_liver_14d_rsem/BM1_liver_14d_rsem.genes.results",
    "/public6/yuanming/paper/publicData/BM2_liver_14d_rsem/BM2_liver_14d_rsem.genes.results",
    "/public6/yuanming/paper/publicData/BM3_liver_14d_rsem/BM3_liver_14d_rsem.genes.results",
    "/public6/yuanming/paper/publicData/BM1_liver_28d_rsem/BM1_liver_28d_rsem.genes.results",
    "/public6/yuanming/paper/publicData/BM2_liver_28d_rsem/BM2_liver_28d_rsem.genes.results",
    "/public6/yuanming/paper/publicData/BM3_liver_28d_rsem/BM3_liver_28d_rsem.genes.results"
)
sampleNames <- c(
    "BM1_0d", "BM2_0d", "BM3_0d",
    "BM1_7d", "BM2_7d", "BM3_7d",
    "BM1_10d", "BM2_10d", "BM3_10d",
    "BM1_14d", "BM2_14d", "BM3_14d",
    "BM1_28d", "BM2_28d", "BM3_28d"
)
condition <- c(
    "0d", "0d", "0d",
    "7d", "7d", "7d",
    "10d", "10d", "10d",
    "14d", "14d", "14d",
    "28d", "28d", "28d"
)
replicate <- c(
    "BM1", "BM2", "BM3",
    "BM1", "BM2", "BM3",
    "BM1", "BM2", "BM3",
    "BM1", "BM2", "BM3",
    "BM1", "BM2", "BM3"
)
comparison <- list(
    c("7d", "0d"),
    c("10d", "0d"), c("10d", "7d"),
    c("14d", "0d"), c("14d", "7d"), c("14d", "10d"),
    c("28d", "0d"), c("28d", "7d"), c("28d", "10d"), c("28d", "14d")
)

count <- readIntoCountMatrix(wholeFilesPath, sampleNames, 5)
col <- generateColDataDF(wholeFilesPath, sampleNames, condition, replicate, ref = rep("mm10", 15))
row <- generateRowDataDF(rownames(count)) %>% tidyr::separate(featureID, into = c("ensembl", "symbol"), sep = "_", remove = F)
vs <- generateComparisonDF(comparison)
tpm <- readIntoCountMatrix(wholeFilesPath, sampleNames, 6)
# save(count, col, row, vs, tpm, file = "/public6/yuanming/paper/Deer/mouse.liver.data.rdata")
# load("/public6/yuanming/paper/Deer/mouse.liver.data.rdata")

mmu <- init_DEanalysis_FromMatrix(
    countMatrix = count, ColDataDF = col, RowDataDF = row, comparisonDF = vs,
    project = "mouse_liver_rnaseq", organism = "mmu", assay.type = "RNA",
    min.count = 15, DE.foldchange = 1, DE.qvalue = 0.05
)
mmu <- appendCount(mmu, countData = tpm, "TPM")
# mmu <- resetAtiveSamples(mmu, drop = c("BM3_7d", "BM2_28d"))
mmu <- diffAnalyze(mmu)

deg.group <- split(deg$entrez, deg$deType.SWA_vs_NC)
input_GSEA_Genes.df <-
    filter(deg, padj.SWA_vs_NC != "NA", symbol != "NA") %>%
    arrange(desc(log2FC.SWA_vs_NC)) %>%
    distinct(symbol, .keep_all = T)
input_GSEA_Genes <- input_GSEA_Genes.df$log2FC.SWA_vs_NC
names(input_GSEA_Genes) <- input_GSEA_Genes.df$symbol
input_GOKEGG_Genes <- list("DEG.down" = unique(na.omit(deg.group$down)), "DEG.up" = unique(na.omit(deg.group$up)))
input_GOKEGG_Keytype <- "ENTREZID"
input_GOKEGG_OrgDb <- org.Mm.eg.db::org.Mm.eg.db
input_KEGG_Org <- "mmu"
input_GSEA_species <- "Mus musculus"


# x <- object
# x[grepl("_brd4", rownames(x), ignore.case = T), ] %>% rowData()
# x[grepl("_brd[234t]$", rownames(x), ignore.case = T), ] %>% assay("TPM")
# x[rowData(x)$symbol %in% c("Brd2", "Brd3", "Brd4"), colData(x)$replicate == "BM1"] %>%
#     rowData() %>%
#     as.data.frame() %>%
#     select(symbol, starts_with("deType"), everything())








### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### mfot
###

wholeFilesPath <- c(
    "/public6/yuanming/paper/publicData/MF1_liver_0d_rsem/MF1_liver_0d_rsem.genes.results",
    "/public6/yuanming/paper/publicData/MF2_liver_0d_rsem/MF2_liver_0d_rsem.genes.results",
    "/public6/yuanming/paper/publicData/MF3_liver_0d_rsem/MF3_liver_0d_rsem.genes.results",
    "/public6/yuanming/paper/publicData/MF1_liver_7d_rsem/MF1_liver_7d_rsem.genes.results",
    "/public6/yuanming/paper/publicData/MF2_liver_7d_rsem/MF2_liver_7d_rsem.genes.results",
    "/public6/yuanming/paper/publicData/MF3_liver_7d_rsem/MF3_liver_7d_rsem.genes.results",
    "/public6/yuanming/paper/publicData/MF1_liver_10d_rsem/MF1_liver_10d_rsem.genes.results",
    "/public6/yuanming/paper/publicData/MF2_liver_10d_rsem/MF2_liver_10d_rsem.genes.results",
    "/public6/yuanming/paper/publicData/MF3_liver_10d_rsem/MF3_liver_10d_rsem.genes.results",
    "/public6/yuanming/paper/publicData/MF1_liver_14d_rsem/MF1_liver_14d_rsem.genes.results",
    "/public6/yuanming/paper/publicData/MF2_liver_14d_rsem/MF2_liver_14d_rsem.genes.results",
    "/public6/yuanming/paper/publicData/MF3_liver_14d_rsem/MF3_liver_14d_rsem.genes.results",
    "/public6/yuanming/paper/publicData/MF1_liver_21d_rsem/MF1_liver_21d_rsem.genes.results",
    "/public6/yuanming/paper/publicData/MF2_liver_21d_rsem/MF2_liver_21d_rsem.genes.results",
    "/public6/yuanming/paper/publicData/MF3_liver_21d_rsem/MF3_liver_21d_rsem.genes.results"
)
sampleNames <- c(
    "MF1_0d", "MF2_0d", "MF3_0d",
    "MF1_7d", "MF2_7d", "MF3_7d",
    "MF1_10d", "MF2_10d", "MF3_10d",
    "MF1_14d", "MF2_14d", "MF3_14d",
    "MF1_28d", "MF2_28d", "MF3_28d"
)
condition <- c(
    "0d", "0d", "0d",
    "7d", "7d", "7d",
    "10d", "10d", "10d",
    "14d", "14d", "14d",
    "28d", "28d", "28d"
)
replicate <- c(
    "BM1", "BM2", "BM3",
    "BM1", "BM2", "BM3",
    "BM1", "BM2", "BM3",
    "BM1", "BM2", "BM3",
    "BM1", "BM2", "BM3"
)
comparison <- list(
    c("7d", "0d"),
    c("10d", "0d"), c("10d", "7d"),
    c("14d", "0d"), c("14d", "7d"), c("14d", "10d"),
    c("28d", "0d"), c("28d", "7d"), c("28d", "10d"), c("28d", "14d")
)
count <- readIntoCountMatrix(wholeFilesPath, sampleNames, 5)
col <- generateColDataDF(wholeFilesPath, sampleNames, condition, replicate, ref = rep("mm10", 15))
row <- generateRowDataDF(rownames(count)) %>% tidyr::separate(featureID, into = c(NULL, "symbol"), sep = "_", remove = F)
vs <- generateComparisonDF(comparison)
tpm <- readIntoCountMatrix(wholeFilesPath, sampleNames, 6)

mfot <- init_DEanalysis_FromMatrix(
    countMatrix = count, ColDataDF = col, RowDataDF = row, comparisonDF = vs,
    project = "mfot_liver_rnaseq", organism = "mmu", assay.type = "RNA",
    min.count = 15, DE.foldchange = 1, DE.qvalue = 0.05
)
mfot <- appendCount(mfot, countData = tpm, "TPM")
mfot <- resetAtiveSamples(mfot, drop = c("MF2_10d", "MF1_14d"))
mfot <- diffAnalyze(mfot)
