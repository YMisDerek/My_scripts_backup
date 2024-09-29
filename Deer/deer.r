rm(list = ls())
suppressMessages(suppressWarnings(library(SummarizedExperiment)))
suppressMessages(suppressWarnings(library(DESeq2)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(methods)))



### ### ### ### ### ### ### ### ### ### ### ### ### ###### ### ### ### ### ### ### ### ### ### ### ### ### ###
### define DEanalysis object
### ### ### ### ### ### ### ### ### ### ### ### ### ###### ### ### ### ### ### ### ### ### ### ### ### ### ###
setClass("diffSet",
    slots = list(
        dds = "DESeqDataSet",
        DE.result = "data.frame",
        DE.stat = "list",
        DE.image = "list"
    )
)

setClass("enrichSet",
    slots = list(
        enrichResults = "list",
        enrich.res = "data.frame",
        enrich.image = "list"
    )
)

setClass("DEanalysis",
    slots = list(
        project = "character",
        organism = "character",
        dataSet = "SummarizedExperiment",
        diffSet = "diffSet",
        enrichSet = "enrichSet"
    ),
    prototype = list(
        project = "NA",
        organism = "NA"
    )
)


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# 定义函数 + 定义方法 + 构建方法
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

check_lengths <- function(...) {
    lengths <- sapply(list(...), length)
    if (length(unique(lengths)) != 1) {
        stop("All input vectors must have the same length.")
    }
}

check_argTypes <- function(args, expected_types) {
    for (i in seq_along(args)) {
        if (!is(args[[i]], expected_types[[i]])) {
            stop(sprintf("Argument '%s' must be of type '%s'.", names(args)[i], expected_types[[i]]))
        }
    }
}


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
readIntoCountMatrix <- function(wholeFilesPath, sampleNames, dataColNum) {
    suppressMessages(suppressWarnings(library(data.table)))
    suppressMessages(suppressWarnings(library(dplyr)))

    check_lengths(wholeFilesPath, sampleNames)

    lapply(seq_along(wholeFilesPath), function(i) {
        fread(wholeFilesPath[i], header = TRUE, select = c(1, dataColNum), col.names = c("featureID", "count"))[, sample := sampleNames[i]]
    }) %>%
        rbindlist() %>%
        dcast(., featureID ~ sample, value.var = "count", fun.aggregate = mean, fill = 0) %>%
        select(featureID, all_of(sampleNames)) %>%
        tibble::column_to_rownames(var = "featureID")
}

generateColDataDF <- function(wholeFilesPath, sampleNames, condition, replicate, ...) {
    check_argTypes(
        args = list(wholeFilesPath, sampleNames, condition, replicate),
        expected_types = c("character", "character", "character", "character")
    )
    check_lengths(wholeFilesPath, sampleNames, condition, replicate)

    coldata <- data.frame(
        sampleLocation = wholeFilesPath,
        sampleName = sampleNames,
        condition = condition,
        replicate = replicate,
        stringsAsFactors = FALSE
    )

    additional_cols <- list(...)
    if (length(additional_cols) > 0) {
        coldata <- cbind(coldata, as.data.frame(additional_cols, stringsAsFactors = FALSE))
    }

    return(coldata)
}

generateRowDataDF <- function(featureID, symbol = NULL, entrez = NULL, ensembl = NULL, ...) {
    fill_missing <- function(vec) {
        if (is.null(vec)) rep(NA, length(featureID)) else vec
    }

    rowdata <- data.frame(
        featureID = featureID,
        symbol = fill_missing(symbol),
        entrez = fill_missing(entrez),
        ensembl = fill_missing(ensembl),
        stringsAsFactors = FALSE
    )

    additional_cols <- list(...)
    if (length(additional_cols) > 0) {
        rowdata <- cbind(rowdata, as.data.frame(additional_cols, stringsAsFactors = FALSE))
    }

    return(rowdata)
}

generateComparisonDF <- function(pairList) {
    if (!is.list(pairList) || !all(sapply(pairList, function(x) length(x) == 2 && all(sapply(x, is.character))))) {
        stop("Input must be a list of pairs: list(c('treat', 'control'))")
    }

    treat <- sapply(pairList, function(x) x[1])
    control <- sapply(pairList, function(x) x[2])

    data.frame(Treat = treat, Control = control, stringsAsFactors = FALSE)
}

initializeDEanalysis_FromMatrix <- function(countMatrix, ColDataDF, RowDataDF, comparisonDF,
    project, organism, assay.type, min.count, DE.foldchange, DE.qvalue) {

    if (assay.type == "RNA") {
        dataSet <- SummarizedExperiment(
            assays = list(count.raw = countMatrix),
            colData = ColDataDF,
            rowData = RowDataDF,
            metadata = list(
                "assay.type" = assay.type,
                "assay.design" = ColDataDF,
                "comparison" = comparisonDF,
                "min.count" = min.count,
                "DE.foldchange" = DE.foldchange,
                "DE.qvalue" = DE.qvalue,
                "low.count.features" = rownames(countMatrix)[rowSums(countMatrix) < min.count],
                "active.samples" = ColDataDF$sampleName
            )
        )
    } else if (assay.type == "ChIP") {
        # rowdatadf要改为granges对象
        dataSet <- RangedSummarizedExperiment(
            assays = list(count.raw = countMatrix),
            colData = ColDataDF,
            rowData = RowDataDF,
            metadata = list(
                "assay.type" = assay.type,
                "assay.design" = ColDataDF,
                "comparison" = comparisonDF,
                "min.count" = min.count,
                "DE.foldchange" = DE.foldchange,
                "DE.qvalue" = DE.qvalue,
                "low.count.features" = rownames(countMatrix)[rowSums(countMatrix) < min.count],
                "active.samples" = ColDataDF$sampleName
            )
        )
    } else {
        stop("Unsupported assay type.")
    }

    DEanalysisObject = new("DEanalysis",
        dataSet = dataSet,
        diffSet = new("diffSet"),
        enrichSet = new("enrichSet")
    )
}

initializeDEanalysis_FromFiles <- function(wholeFilesPath, sampleNames, dataColNum, condition, replicate, pairList,
    project, organism, assay.type, min.count, DE.foldchange, DE.qvalue) {
    
    countMatrix <- readIntoCountMatrix(wholeFilesPath, sampleNames, dataColNum)
    ColDataDF <- generateColDataDF(wholeFilesPath, sampleNames, condition, replicate)
    RowDataDF <- generateRowDataDF(rownames(countMatrix))
    comparisonDF <- generateComparisonDF(pairList)

    if (assay.type == "RNA") {
        dataSet <- SummarizedExperiment(
            assays = list(count.raw = countMatrix),
            colData = ColDataDF,
            rowData = RowDataDF,
            metadata = list(
                "assay.type" = assay.type,
                "assay.design" = ColDataDF,
                "comparison" = comparisonDF,
                "min.count" = min.count,
                "DE.foldchange" = DE.foldchange,
                "DE.qvalue" = DE.qvalue,
                "low.count.features" = rownames(countMatrix)[rowSums(countMatrix) < min.count],
                "active.samples" = ColDataDF$sampleName
            )
        )
    } else if (assay.type == "ChIP") {
        # rowdatadf要改为granges对象
        dataSet <- RangedSummarizedExperiment(
            assays = list(count.raw = countMatrix),
            colData = ColDataDF,
            rowData = RowDataDF,
            metadata = list(
                "assay.type" = assay.type,
                "assay.design" = ColDataDF,
                "comparison" = comparisonDF,
                "min.count" = min.count,
                "DE.foldchange" = DE.foldchange,
                "DE.qvalue" = DE.qvalue,
                "low.count.features" = rownames(countMatrix)[rowSums(countMatrix) < min.count],
                "active.samples" = ColDataDF$sampleName
            )
        )
    } else {
        stop("Unsupported assay type.")
    }

    DEanalysisObject = new("DEanalysis",
        dataSet = dataSet,
        diffSet = new("diffSet"),
        enrichSet = new("enrichSet")
    )
}

add_countData <- function(DEanalysisObject, data, dataType) {
    if (dataType == "TPM") {
        assays(DEanalysisObject@dataSet)$TPM <- data
    } else if (dataType == "FPKM") {
        assays(DEanalysisObject@dataSet)$FPKM <- data
    } else {
        stop("only 'TPM' or 'FPKM' recepted...")
    }
}

refilterDataSet <- function(DEanalysisObject, min.count) {
    metadata(DEanalysisObject@dataSet)$min.count <- min.count
    metadata(DEanalysisObject@dataSet)$low.count.features <-
        rownames(DEanalysisObject@dataSet)[
            rowSums(assay(DEanalysisObject@dataSet, "count.raw")) < metadata(DEanalysisObject@dataSet)$min.count
        ]
}

# DiffAnalyze <- function(DEanalysisObject) {
#     suppressWarnings(suppressMessages(library(DESeq2)))
#     suppressWarnings(suppressMessages(library(dplyr)))
#     suppressWarnings(suppressMessages(library(tidyr)))
#     suppressWarnings(suppressMessages(library(tibble)))
#     suppressWarnings(suppressMessages(library(purrr)))
#     suppressWarnings(suppressMessages(library(stringr)))
#     suppressWarnings(suppressMessages(library(ggplot2)))
#     suppressWarnings(suppressMessages(library(lay)))
#     suppressWarnings(suppressMessages(library(ggrastr)))

#     filter <- metadata(DEanalysisObject)$min.count
#     cutoff <- metadata(DEanalysisObject)$DE.foldchange
#     signif <- metadata(DEanalysisObject)$DE.qvalue
#     projectName <- DEanalysisObject@project

#     # message("========================diff analysis: DEseq2==========================")
#     # message("Working Dir": getwd())
#     # message("filter out counts less than: ", filter)
#     # message("significant cutoff         : ", signif)
#     # message("log2foldchange cutoff      : ", cutoff)
#     # message("result: ", "\t", paste0(projectName, "*"))
#     # message("=======================================================================\n")

#     # DEobj@diff.Object@dds <- DESeqDataSetFromMatrix(
#     #     countData = DEobj@dataSet@count.raw,
#     #     colData = DEobj@dataSet@colData,
#     #     design = ~condition
#     # )
#     # 过滤低count
#     dds <- DESeq(dds, fitType = "parametric")
#     rld <- rlog(dds, blind = FALSE)
#     counts.rlog <- assay(rld) %>% as_tibble(., rownames = "ID")
# }

# EnrichAnalyze <- function(DEanalysisObject) {}







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
    "/public6/yuanming/paper/publicData/BM1_liver_7d_rsem/BM1_liver_7d_rsem.genes.results"
)
sampleNames <- c("BM1_0d", "BM1_7d")
condition <- c("0d", "7d")
replicate <- c("1", "1")
comparison <- list(c("7d", "0d"))

count <- readIntoCountMatrix(wholeFilesPath, sampleNames, 5)
tpm <- readIntoCountMatrix(wholeFilesPath, sampleNames, 6)
col <- generateColDataDF(wholeFilesPath, sampleNames, condition, replicate, ref = c("mm10", "mm10"))
row <- generateRowDataDF(rownames(count)) %>% tidyr::separate(featureID, into = c("ensembl", "symbol"), sep = "_", remove = F)
vs = generateComparisonDF(comparison)

x = initializeDEanalysis_FromMatrix(
    countMatrix = count, ColDataDF = col, RowDataDF = row, comparisonDF = vs,
    project = 'mytest', organism = 'mouse', assay.type = 'RNA', 
    min.count = 10, DE.foldchange = 1, DE.qvalue = 0.05
)
