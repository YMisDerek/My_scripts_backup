### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### class
###

setClass("DEanalysis",
    slots = list(
        project = "character_OR_NULL",
        organism = "character_OR_NULL",
        dataSet = "RangedSummarizedExperiment",
        diffSet = "diffSet",
        enrichSet = "enrichSet"
    ),
    prototype = list(
        project = "NA",
        organism = "NA"
    )
)


initializeDEanalysis_FromMatrix <- function(countMatrix, ColDataDF, RowDataDF, comparisonDF,
                                            project, organism = c("hsa", "mmu", "mfot"),
                                            assay.type = c("RNA", "ChIP"),
                                            min.count = 0, DE.foldchange = 1, DE.qvalue = 0.05) {
    assay.type <- match.arg(assay.type, choices = c("RNA", "ChIP"))
    organism <- match.arg(organism, choices = c("hsa", "mmu", "mfot"))

    dataSet <- createDataSet(
        countMatrix, ColDataDF, RowDataDF, comparisonDF,
        assay.type, assay.design, min.count, DE.foldchange, DE.qvalue
    )

    new("DEanalysis",
        project = project,
        organism = organism,
        dataSet = dataSet,
        diffSet = new("diffSet"),
        enrichSet = new("enrichSet")
    )
}



### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### manipulate
###

setGeneric("appendCount", signature = "obj", function(obj, countData, dataName, ...) standardGeneric("appendCount"))

.append_countData_to_dataSet <- function(obj, countData, dataName = c("TPM", "FPKM", "rlog")) {
    dataName <- match.arg(dataName, choices = c("TPM", "FPKM", "rlog"))
    assays(obj@dataSet)[[dataName]] <- countData
    return(obj)
}

setMethod("appendCount", "DEanalysis", .append_countData_to_dataSet)


refilterDataSet <- function(obj, min.count) {
    metadata(obj@dataSet)$min.count <- min.count
    metadata(obj@dataSet)$actifeave.features <-
        rownames(obj@dataSet)[
            rowSums(assay(obj@dataSet, "count.raw")) > metadata(obj@dataSet)$min.count
        ]
    return(obj)
}

reActiveSample <- function(obj, deactivate.sample) {
    all.sample <- colnames(obj@dataSet)
    keep <- !(all.sample %in% deactivate.sample)
    obj@dataSet@metadata$active.samples <- all.sample[keep]
    obj@dataSet@metadata$de.active.samples <- as.vector(deactivate.sample)
    return(obj)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### show
###

setMethod(
    "show", "DEanalysis",
    function(object) {
        cat("class:", class(object), "\n")
        cat("project:", object@project, "\n")
        cat("organism:", object@organism, "\n")

        ### dataSet
        cat("\n.. @dataSet\n")
        coolcat("assays(%d): %s\n", assayNames(obj@dataSet))
        cat("dim:", dim(object@dataSet), "\n")
        coolcat("acitve features(%d): %s\n", metadata(object@dataSet)$active.features)
        coolcat("acitve samples(%d): %s\n", metadata(object@dataSet)$active.samples)
        coolcat("rowData names(%d): %s\n", names(rowData(object@dataSet)))
        coolcat("colData names(%d): %s\n", names(colData(object@dataSet)))
        cat(
            "current threshold:",
            paste0(
                "min.count(", metadata(object@dataSet)$min.count, "), ",
                "DE.foldchange(", metadata(object@dataSet)$DE.foldchange, "), ",
                "DE.qvalue(", metadata(object@dataSet)$DE.qvalue, ")\n"
            )
        )

        ### diffSet
        cat("\n.. @diffSet\n")
    }
)
