### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### diffAnalyze flow
###

diffAnalyze <- function(object, clean = TRUE) {
    if (!is(object, "DEanalysis")) stop("input object not DEanalysis")


    ### get input value
    myObj <- object[
        rownames(object) %in% metadata(object)$active.features,
        object$sampleName %in% metadata(object)$active.samples
    ]
    outBaseName <- metadata(myObj)$project
    comparisonDF <- metadata(myObj)$comparisonDF
    signif <- metadata(myObj)$DE.qvalue
    cutoff <- metadata(myObj)$DE.foldchange


    ### deseq
    ### save diffSet@DE.dds
    stopifnot("counts" %in% names(assays(myObj)))
    stopifnot("condition" %in% colnames(colData(myObj)))
    dds <- DESeqDataSetFromMatrix(countData = round(assay(myObj, "counts")), colData = colData(myObj), design = ~condition)
    dds <- DESeq(dds, fitType = "parametric")
    object@diffSet@DE.dds <- dds


    ### save diff result into rowData
    DEresults <- .get_DESeq_result_all(dds, comparisonDF, signif, cutoff, clean = TRUE)
    rowData(object) <- left_join(as.data.frame(rowData(object)), DEresults, by = "featureID")


    ### save rlog into assays
    rld <- rlog(dds, blind = FALSE)
    counts.rlog <- assay(rld) %>% as_tibble(., rownames = "ID")

    rlogData <- as.data.frame(matrix(NA, nrow = nrow(object), ncol = ncol(object), dimnames = dimnames(object)))
    rows.idx <- match(rownames(myObj), rownames(object))
    cols.idx <- match(colnames(myObj), colnames(object))
    rlogData[rows.idx, cols.idx] <- counts.rlog[, -1]

    object <- appendCount(object, countData = rlogData, "rlog")


    ### save diffSet@DE.image
    DE.image <- .QC.assay(rld, counts.rlog, myObj$sampleName)
    object@diffSet@DE.image <- DE.image


    ### save diffSet@DE.stat
    DE.num <- map_df(1:nrow(comparisonDF), function(i) {
        treat <- comparisonDF$Treat[i]
        control <- comparisonDF$Control[i]
        stat <- table(DEresults[, paste0("deType", ".", treat, "_vs_", control)]) %>%
            as.list() %>%
            as.data.frame()
        stat$comparison <- paste0(treat, "_vs_", control)
        return(stat)
    })
    object@diffSet@DE.stat <- DE.num[, c("comparison", "Down", "None", "Up")]


    return(object)
}

.QC.assay <- function(rld, counts.rlog, sampleName, ...) {
    fig.pca <-
        plotPCA(rld) +
        background +
        theme(legend.position = "none", aspect.ratio = 1) +
        coord_cartesian() +
        ggrepel::geom_label_repel(aes(label = name), size = 8 / .pt)
    fig.count <-
        pivot_longer(counts.rlog, -ID) %>%
        ggplot(aes(x = factor(name, levels = sampleName), y = value)) +
        geom_boxplot(notch = T, fill = "grey") +
        labs(x = NULL, y = "normalized signal") +
        background +
        theme(axis.text.x = element_text(angle = 90, hjust = 1), aspect.ratio = 1)
    samplesCorrDF <-
        Hmisc::rcorr(as.matrix(counts.rlog[, -1]))$r %>%
        as.data.frame() %>%
        rownames_to_column() %>%
        pivot_longer(cols = -1) %>%
        rename_with(~ c("name1", "name2", "rho"))
    fig.corr <-
        ggplot(samplesCorrDF) +
        geom_tile(aes(factor(name1, levels = sampleName),
            factor(name2, levels = sampleName),
            fill = rho
        )) +
        geom_text(aes(factor(name1, levels = sampleName),
            factor(name2, levels = sampleName),
            label = round(rho, 2)
        ), size = 2) +
        theme_classic() +
        theme(
            axis.text.x = element_text(angle = 30, hjust = 1),
            axis.text = element_text(color = "black"),
            axis.title = element_blank(),
            aspect.ratio = 1
        ) +
        scale_fill_gradient(low = "white", high = "#1b98e0")

    return(list("PCA" = fig.pca, "normBar" = fig.count, "sampleCorr" = fig.corr))
}

.get_DESeq_result_all <- function(dds, comparisonDF, signif, cutoff, clean) {
    results_list <- lapply(1:nrow(comparisonDF), function(i) {
        treat <- comparisonDF$Treat[i]
        control <- comparisonDF$Control[i]
        .get_DESeq_result(dds, treat, control, signif, cutoff, clean)
    })
    reduce(results_list, full_join, by = "featureID")
}

.get_DESeq_result <- function(dds, treat, control, signif, cutoff, clean) {
    stopifnot(is.logical(clean))

    message("extracting treat/control.. \t", treat, " vs ", control)
    DEresult <-
        results(dds, contrast = c("condition", treat, control)) %>%
        as_tibble(rownames = "featureID") %>%
        mutate(deType = case_when(
            padj < signif & log2FoldChange > cutoff ~ "Up",
            padj < signif & log2FoldChange < cutoff * (-1) ~ "Down",
            TRUE ~ "None"
        ))
    colnames(DEresult) <- c("featureID", paste0(colnames(DEresult)[-1], paste0(".", treat, "_vs_", control)))

    if (clean) {
        DEresult <- DEresult[, -c(2, 4, 5)]
    }
    DEresult
}

# .stat_plot <- function(DEresult, treat, control, counts.rlog, ...) {
#     changeNum <-
#         dplyr::count(DEresult, deType, name = "Freq") %>%
#         mutate(tag = paste0(deType, ": ", Freq))
#     xlab_range <- DEresult$log2FC %>%
#         abs() %>%
#         max() %>%
#         ceiling()
#     ylab_range <- DEresult$padj %>%
#         -log10(.) %>%
#         max() %>%
#         ceiling()
#     volc <-
#         ggplot(DEresult) +
#         geom_point_rast(aes(x = log2FC, y = -log10(padj), color = deType), size = 0.3) +
#         geom_vline(xintercept = c(0), color = "#990000", linetype = "dashed") +
#         scale_color_manual(
#             values = c("Down" = "#4d97cd", "None" = "#dee2e6", "Up" = "#c74546"),
#             labels = setNames(changeNum$tag, changeNum$deType),
#             name = paste0(treat, "\nvs\n", control)
#         ) +
#         scale_x_continuous(
#             limits = c(-xlab_range, xlab_range),
#             breaks = seq(-xlab_range, xlab_range, by = 1)
#         ) +
#         background +
#         theme(aspect.ratio = 1) +
#         guides(color = guide_legend(override.aes = list(size = 5)))
#     sca <-
#         dplyr::mutate(counts.rlog,
#             treat_norm = lay(across(contains(treat)), mean),
#             control_norm = lay(across(contains(control)), mean)
#         ) %>%
#         dplyr::select(ID, treat_norm, control_norm) %>%
#         dplyr::left_join(., dplyr::select(DEresult, ID, deType), by = "ID") %>%
#         ggplot(aes(x = control_norm, y = treat_norm, color = deType)) +
#         geom_point_rast(size = 0.2, alpha = 0.7) +
#         geom_abline(slope = 1, intercept = 0, lty = 2) +
#         labs(x = paste0("Intensity of ", control), y = paste0("Intensity of ", treat), color = "deType") +
#         scale_color_manual(values = c("Down" = "#4d97cd", "None" = "#dee2e6", "Up" = "#c74546")) +
#         scale_x_continuous(breaks = c(0, 3, 6, 9, 12, 15), limits = c(0, 15)) +
#         scale_y_continuous(breaks = c(0, 3, 6, 9, 12, 15), limits = c(0, 15)) +
#         background +
#         guides(color = guide_legend(override.aes = list(size = 5))) +
#         theme(legend.position = "none", aspect.ratio = 1)
#     print(ggpubr::ggarrange(volc, sca, nrow = 1, ncol = 2, align = "h", common.legend = T, legend = "right"))
# }

# plot.volc = function(treat, control)
