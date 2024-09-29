rm(list = ls())
suppressMessages(suppressWarnings(library(devtools)))

setwd("D:\\300_学习\\340_博士后\\Deer\\deer")
load_all()

load("../test.rdata")

obj <- initializeDEanalysis_FromMatrix(
    countMatrix = count, ColDataDF = col, RowDataDF = row, comparisonDF = vs,
    project = "mytest", organism = "mmu", assay.type = "RNA",
    min.count = 40, DE.foldchange = 1, DE.qvalue = 0.05
)
obj = appendCount(obj, countData = tpm, 'TPM')
obj



# class(obj@dataSet)
# x = obj@dataSet
# x[grepl("Gnai3", rownames(x)), ] %>% rowData()
# x[rowData(x)$symbol == "Gnai3", 'BM3_7d'] %>% assay()
