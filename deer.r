rm(list = ls())
suppressMessages(suppressWarnings(library(SummarizedExperiment)))
suppressMessages(suppressWarnings(library(DESeq2)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(methods)))



### ### ### ### ### ### ### ### ### ### ### ### ### ###### ### ### ### ### ### ### ### ### ### ### ### ### ###### ### ### ### ### ### ### ### ### ### ### ### ### ###
### define DEanalysis object
### ### ### ### ### ### ### ### ### ### ### ### ### ###### ### ### ### ### ### ### ### ### ### ### ### ### ###### ### ### ### ### ### ### ### ### ### ### ### ### ###

# 基类BaseData
setClass("BaseData",
    slots = list(
        data = "ANY",
        low.count.features = "character",
        active.samples = "character"
    ),
    prototype = list(
        low.count.features = "",
        active.samples = ""
    )
)

setClass("RNA.data",
    contains = "BaseData",
    slots = list(data = "SummarizedExperiment")
)

setClass("ChIP.data",
    contains = "BaseData",
    slots = list(data = "RangedSummarizedExperiment")
)

setClass("diff.Object",
    slots = list(
        dds = "DESeqDataSet",
        DE.result = "data.frame",
        DE.stat = "list",
        DE.image = "list"
    )
)

setClass("enrich.Object",
    slots = list(
        enrichResults = "list",
        enrich.res = "data.frame",
        enrich.image = "list"
    )
)

setClass("assayInfo",
    slots = list(
        project.name = "character",
        organism = "character",
        assay.type = "character",
        min.count = "numeric",
        DE.foldchange = "numeric",
        DE.qvalue = "numeric",
        assay.design = "data.frame",
        comparison = "data.frame"
    ), prototype = list(
        min.count = 0,
        DE.foldchange = 1,
        DE.qvalue = 0.05
    )
)

setClass("DEanalysis",
    slots = list(
        assayInfo = "assayInfo",
        data.Object = "BaseData",
        diff.Object = "diff.Object",
        enrich.Object = "enrich.Object"
    )
)


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# 定义函数 + 定义方法 + 构建方法
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

# 读取rsem定量结果 + 合并数据
readIntoCountMatrix <- function(wholeFilesPath, sampleNames, countType = "count", countColNum = NULL) {
    suppressMessages(suppressWarnings(library(data.table)))
    suppressMessages(suppressWarnings(library(dplyr)))

    read_counts <- function(file_path, countType) {
        offset <- if (!is.null(countColNum)) {
            countColNum
        } else {
            switch(countType, "TPM" = 5, "FPKM" = 6, "count" = 4, stop('Invalid countType. Must be "TPM", "FPKM", or "count".'))
        }
        fread(file_path, header = TRUE, select = c(1, offset + 1), col.names = c("featureID", "count"))
    }

    create_data_matrix <- function(wholeFilesPath, sampleNames, countType) {
        lapply(seq_along(wholeFilesPath), function(i) {
            read_counts(wholeFilesPath[i], countType)[, sample := sampleNames[i]]
        }) %>%
            rbindlist() %>%
            dcast(., featureID ~ sample, value.var = "count", fun.aggregate = mean, fill = 0) %>%
            select(featureID, all_of(sampleNames)) %>%
            tibble::column_to_rownames(var = "featureID")
    }

    create_data_matrix(wholeFilesPath, sampleNames, countType)
}

# 生成coldata
generateColDataTable <- function(wholeFilesPath = "", sampleNames = "", condition = "", replicate = "", ...) {
    check_argTypes <- function(args, expected_types) {
        for (i in seq_along(args)) {
            if (!is(args[[i]], expected_types[[i]])) {
                stop(sprintf("Argument '%s' must be of type '%s'.", names(args)[i], expected_types[[i]]))
            }
        }
    }
    check_argTypes(
        args = list(wholeFilesPath = wholeFilesPath, sampleNames = sampleNames, condition = condition, replicate = replicate),
        expected_types = c("character", "character", "character", "character")
    )

    check_lengths <- function(...) {
        lengths <- sapply(list(...), length)
        if (length(unique(lengths)) != 1) {
            stop("All input vectors must have the same length.")
        }
    }
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
        additional_df <- as.data.frame(additional_cols, stringsAsFactors = FALSE)
        coldata <- cbind(coldata, additional_df)
    }

    return(coldata)
}

# 生成rowdata
generateRowDataTable <- function(featureID, symbol = NULL, entrez = NULL, ensembl = NULL, ...) {
    fill_missing <- function(vec) { if (is.null(vec)) rep(NA, length(featureID)) else vec }

    rowdata <- data.frame(
        featureID = featureID,
        symbol = fill_missing(symbol),
        entrez = fill_missing(entrez),
        ensembl = fill_missing(ensembl),
        stringsAsFactors = FALSE
    )

    additional_cols <- list(...)
    if (length(additional_cols) > 0) {
        additional_df <- as.data.frame(additional_cols, stringsAsFactors = FALSE)
        rowdata <- cbind(rowdata, additional_df)
    }

    return(rowdata)
}

# 生成comparison
generateComparisonTable <- function(pairList) {
    if (!is.list(pairList) || !all(sapply(pairList, function(x) length(x) == 2 && all(sapply(x, is.character))))) {
        stop("Input must be a list of pairs: list(c('treat', 'control'))")
    }

    treat <- sapply(pairList, function(x) x[1])
    control <- sapply(pairList, function(x) x[2])

    data.frame(Treat = treat, Control = control, stringsAsFactors = FALSE)
}

# 创建数据对象
createDataObject <- function(assay.type, count.matrix, ColDataTable, RowDataTable) {
    if (assay.type == "RNA") {
        data <- new("RNA.data",
            data = SummarizedExperiment(
                assays = list(count.raw = count.matrix),
                colData = ColDataTable,
                rowData = RowDataTable,
                metadata = list("current.min.count" = NULL, "featureID.type" = NULL)
            )
        )
    } else if (assay.type == "ChIP") {
        data <- new("ChIP.data",
            data = RangedSummarizedExperiment(
                assays = list(count.raw = count.matrix, TPM = NULL, FPKM = NULL),
                colData = ColDataTable,
                rowData = RowDataTable,
                metadata = list("current.min.count" = NULL, "featureID.type" = NULL)
            )
        )
    } else {
        stop("Unsupported data type.")
    }
    return(data)
}

# 创建DEanalysis对象
creatDEanalysisObject <- function(project, organism,
                                  count.raw, ColDataTable, RowDataTable, assay.type = "RNA",
                                  compareList,
                                  min.count = 0,
                                  DE.foldchange = 1, DE.qvalue = 0.05) {
    assayInfo <- new("assayInfo",
        project.name = project,
        organism = organism,
        assay.type = assay.type,
        min.count = min.count,
        DE.foldchange = DE.foldchange,
        DE.qvalue = DE.qvalue,
        assay.design = ColDataTable,
        comparison = generateComparisonTable(compareList)
    )
    data.Object <- createDataObject(assay.type, count.raw, ColDataTable, RowDataTable)

    primary_DEanalysis <- new("DEanalysis",
        assayInfo = assayInfo,
        data.Object = data.Object,
        diff.Object = new("diff.Object"),
        enrich.Object = new("enrich.Object")
    )

    return(primary_DEanalysis)
}

add_countData <- function(DEanalysisObject, data, dataType) {
    if (dataType == "TPM") {
        assays(DEanalysisObject@data.Object@data)$TPM <- data
    } else if (dataType == "FPKM") {
        assays(DEanalysisObject@data.Object@data)$FPKM <- data
    } else {
        stop("only 'TPM' or 'FPKM' recepted...")
    }
}

DiffAnalyze <- function(DEanalysisObject) {
    suppressWarnings(suppressMessages(library(DESeq2)))
    suppressWarnings(suppressMessages(library(dplyr)))
    suppressWarnings(suppressMessages(library(tidyr)))
    suppressWarnings(suppressMessages(library(tibble)))
    suppressWarnings(suppressMessages(library(purrr)))
    suppressWarnings(suppressMessages(library(stringr)))
    suppressWarnings(suppressMessages(library(ggplot2)))
    suppressWarnings(suppressMessages(library(lay)))
    suppressWarnings(suppressMessages(library(ggrastr)))

    DEobj <- DEanalysisObject
    filter <- DEobj@assayInfo@min.count
    cutoff <- DEobj@assayInfo@DE.foldchange
    signif <- DEobj@assayInfo@DE.qvalue
    projecName <- DEobj@assayInfo@project.name

    # message("========================diff analysis: DEseq2==========================")
    # message("Working Dir": getwd())
    # message("filter out counts less than: ", filter)
    # message("significant cutoff         : ", signif)
    # message("log2foldchange cutoff      : ", cutoff)
    # message("result: ", "\t", paste0(outputFileName, "*"))
    # message("=======================================================================\n")

    # dds <- DESeqDataSetFromMatrix(countData = counts, colData = colData, design = ~condition)
    dds <- DESeqDataSet(countData = counts, colData = colData, design = ~condition)
    dds <- DESeq(dds, fitType = "parametric")
    rld <- rlog(dds, blind = FALSE)
    counts.rlog <- assay(rld) %>% as_tibble(., rownames = "ID")
}

EnrichAnalyze <- function(DEanalysisObject) {}







### 定义方法
# 输入基因list，绘制表达量TPM柱状图 附加样品散点
# 根据PCA结果，修改active.sample进行后续分析

# ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
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

count <- readIntoCountMatrix(wholeFilesPath, sampleNames, "count")
tpm <- readIntoCountMatrix(wholeFilesPath, sampleNames, "TPM")
col <- generateColDataTable(wholeFilesPath, sampleNames, condition, replicate, ref = c("mm10", "mm10"))
row <- generateRowDataTable(rownames(count)) %>% tidyr::separate(featureID, into = c("ensembl", "symbol"), sep = "_", remove = F)

x = creatDEanalysisObject(
    project = "test", organism = "mouse",
    count.raw = count, ColDataTable = col, RowDataTable = row,
    compareList = comparison, min.count = 50
)
