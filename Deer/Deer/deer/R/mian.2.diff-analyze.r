DiffAnalyze <- function(obj) {
    cutoff <- metadata(obj@dataSet)$DE.foldchange
    signif <- metadata(obj@dataSet)$DE.qvalue
    comparison <- metadata(obj@dataSet)$comparison
    projectName <- obj@project

    message("\n========================diff analysis: DEseq2==========================")
    message(paste0("Working Dir                : ", getwd()))
    message("filter out counts less than: ", metadata(obj@dataSet)$min.count)
    message("significant cutoff         : ", signif)
    message("log2foldchange cutoff      : ", cutoff)
    message("result                     : ", paste0(projectName, "*"))
    message("=======================================================================\n")

    dds <- DESeqDataSetFromMatrix(
        countData = round(assay(obj@dataSet, "count.raw")),
        colData = colData(obj@dataSet),
        design = ~condition
    )
    dds <- dds[
        metadata(obj@dataSet)$active.features,
        metadata(obj@dataSet)$active.samples
    ]

    dds <- DESeq(dds, fitType = "parametric")
    rld <- rlog(dds, blind = FALSE)
    counts.rlog <- assay(rld) %>%
        as_tibble(., rownames = "featureID") %>%
        left_join(
            rowData(obj@dataSet)[, "featureID", drop = F] %>% as.data.frame(),
            .,
            by = "featureID"
        ) %>%
        column_to_rownames(var = "featureID")
    tmp <- data.frame(matrix(nrow = nrow(obj@dataSet), ncol = length(obj@dataSet@metadata$de.active.samples)))
    colnames(tmp) <- obj@dataSet@metadata$de.active.samples
    counts.rlog <- cbind(counts.rlog, tmp) %>% select(all_of(colnames(obj@dataSet)))

    sapply(seq_len(nrow(comparison)), function(i) {
        treat <- comparison$Treat[i]
    })

    obj@diffSet@dds <- dds
    assays(obj@dataSet)$count.rlog <- counts.rlog
}