### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### 定义方法
# 输入基因list，绘制表达量TPM柱状图 附加样品散点
# 根据PCA结果，修改active.sample进行后续分析

# 输出xls文件
# featureID，symbol，entrez，ensembl，log2FC.DE1, qvalue.DE1, DE2, ...
# featureID, TPM, count.norm, FPKM, count.raw
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
wholeFilesPath <- c(
    "/public6/yuanming/paper/publicData/BM1_liver_0d_rsem/BM1_liver_0d_rsem.genes.results",
    "/public6/yuanming/paper/publicData/BM2_liver_0d_rsem/BM2_liver_0d_rsem.genes.results",
    "/public6/yuanming/paper/publicData/BM3_liver_0d_rsem/BM3_liver_0d_rsem.genes.results",
    "/public6/yuanming/paper/publicData/BM1_liver_7d_rsem/BM1_liver_7d_rsem.genes.results",
    "/public6/yuanming/paper/publicData/BM2_liver_7d_rsem/BM2_liver_7d_rsem.genes.results",
    "/public6/yuanming/paper/publicData/BM3_liver_7d_rsem/BM3_liver_7d_rsem.genes.results"
)
sampleNames <- c("BM1_0d", "BM2_0d", "BM3_0d", "BM1_7d", "BM2_7d", "BM3_7d")
condition <- c("0d", "0d", "0d", "7d", "7d", "7d")
replicate <- c("BM1", "BM2", "BM3", "BM1", "BM2", "BM3")
comparison <- list(c("7d", "0d"))

count <- readIntoCountMatrix(wholeFilesPath, sampleNames, 5)
tpm <- readIntoCountMatrix(wholeFilesPath, sampleNames, 6)
col <- generateColDataDF(wholeFilesPath, sampleNames, condition, replicate, ref = c("mm10", "mm10"))
row <- generateRowDataDF(rownames(count)) %>% tidyr::separate(featureID, into = c("ensembl", "symbol"), sep = "_", remove = F)
vs <- generateComparisonDF(comparison)

obj <- initializeDEanalysis_FromMatrix(
    countMatrix = count, ColDataDF = col, RowDataDF = row, comparisonDF = vs,
    project = "mytest", organism = "mmu", assay.type = "RNA",
    min.count = 40, DE.foldchange = 1, DE.qvalue = 0.05
)
