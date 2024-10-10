### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### DEanalysis object
###

setClass("DEanalysis",
    contains = "RangedSummarizedExperiment",
    slots = list(
        diffSet = "diffSet",
        enrichSet = "enrichSet"
    )
)

setValidity("DEanalysis", function(object) {
    metaList <- list(
        "project",
        "organism",
        "assay.type",
        "comparisonDF",
        "min.count",
        "DE.foldchange",
        "DE.qvalue",
        "active.features",
        "active.samples"
    )

    metadata_names <- names(metadata(object))

    missing_vars <- setdiff(metaList, metadata_names)
    if (length(missing_vars) > 0) {
        return(paste("缺少变量:", paste(missing_vars, collapse = ", ")))
    }

    extra_vars <- setdiff(metadata_names, metaList)
    if (length(extra_vars) > 0) {
        return(paste("多余的变量:", paste(extra_vars, collapse = ", ")))
    }

    return(TRUE)
})


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### constructor
###

init_DEanalysis <- function(rse) {
    new("DEanalysis", rse, diffSet = new("diffSet"), enrichSet = new("enrichSet"))
}

init_DEanalysis_FromMatrix <- function(countMatrix, ColDataDF, RowDataDF, comparisonDF,
                                       project, organism = c("hsa", "mmu", "mfot"),
                                       assay.type = c("RNA", "ChIP"),
                                       min.count = 0, DE.foldchange = 1, DE.qvalue = 0.05) {
    assay.type <- match.arg(assay.type, choices = c("RNA", "ChIP"))
    organism <- match.arg(organism, choices = c("hsa", "mmu", "mfot"))

    rse <- SummarizedExperiment(
        assays = list(counts = countMatrix),
        colData = ColDataDF,
        rowData = RowDataDF,
        metadata = list(
            "project" = project,
            "organism" = organism,
            "assay.type" = assay.type,
            "comparisonDF" = comparisonDF,
            "min.count" = min.count,
            "DE.foldchange" = DE.foldchange,
            "DE.qvalue" = DE.qvalue,
            "active.features" = rownames(countMatrix)[rowSums(countMatrix) > min.count],
            "active.samples" = ColDataDF$sampleName
        )
    )
    rse <- as(rse, "RangedSummarizedExperiment")

    init_DEanalysis(rse)
}








### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### method
###

# append matrix to assay
setGeneric("appendCount", signature = "obj", function(obj, countData, dataName, ...) standardGeneric("appendCount"))

.append_countData <- function(obj, countData, dataName = c("TPM", "FPKM", "rlog")) {
    dataName <- match.arg(dataName, choices = c("TPM", "FPKM", "rlog"))
    assays(obj)[[dataName]] <- countData
    return(obj)
}

setMethod("appendCount", "DEanalysis", .append_countData)

# reset filtering
resetFilter <- function(obj, min.count) {
    stopifnot(is(obj, "DEanalysis"))

    metadata(obj)$min.count <- min.count
    metadata(obj)$active.features <- rownames(obj)[rowSums(assay(obj, "counts")) >= metadata(obj)$min.count]
    return(obj)
}

# Reset active samples by removing specified samples
# 还需要补充更新comparison列表，可以配合condition检验
resetAtiveSamples <- function(obj, drop) {
    stopifnot(is(obj, "DEanalysis"))
    stopifnot(is.character(drop))

    all.samples <- colnames(obj)
    stopifnot(all(drop %in% all.samples))
    keep <- !(all.samples %in% drop)

    metadata(obj)$active.samples <- all.samples[keep]
    metadata(obj)$excluded.samples <- drop

    return(obj)
}



### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### show
###

setMethod("show", "DEanalysis", function(object) {
    cat("class:", class(object), "\n")
    cat("project:", metadata(object)$project, "\n")
    cat("metadata:", names(metadata(object)), "\n")

    coolcat("assays(%d): %s\n", assayNames(object))
    cat("dim:", dim(object), "\n")
    coolcat("acitve features(%d): %s\n", metadata(object)$active.features)
    coolcat("acitve samples(%d): %s\n", metadata(object)$active.samples)

    coolcat("rowData names(%d): %s\n", names(rowData(object)))
    coolcat("colData names(%d): %s\n", names(colData(object)))
    cat(
        "current threshold:",
        paste0(
            "min.count(", metadata(object)$min.count, "), ",
            "DE.foldchange(", metadata(object)$DE.foldchange, "), ",
            "DE.qvalue(", metadata(object)$DE.qvalue, ")\n"
        )
    )

    ### diffSet
    print(object@diffSet)

    ### enrichSet
    # object@enrichSet
})
