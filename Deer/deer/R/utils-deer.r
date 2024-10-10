### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### check
###

.check_lengths <- function(...) {
    lengths <- sapply(list(...), length)
    if (length(unique(lengths)) != 1) {
        stop("All input vectors must have the same length.")
    }
}

.check_argTypes <- function(args, expected_types) {
    for (i in seq_along(args)) {
        if (!is(args[[i]], expected_types[[i]])) {
            stop(sprintf("Argument '%s' must be of type '%s'.", names(args)[i], expected_types[[i]]))
        }
    }
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### main
###

readIntoCountMatrix <- function(wholeFilesPath, sampleNames, dataColNum) {
    .check_lengths(wholeFilesPath, sampleNames)
    .check_argTypes(
        list(wholeFilesPath, sampleNames, dataColNum),
        c("character", "character", "numeric")
    )

    lapply(seq_along(wholeFilesPath), function(i) {
        # fread(wholeFilesPath[i], header = TRUE, select = c(1, dataColNum), col.names = c("featureID", "count"))[, sample := sampleNames[i]]
        fread(wholeFilesPath[i], header = TRUE, select = c(1, dataColNum), col.names = c("featureID", "count")) %>%
            mutate(sample = sampleNames[i])
    }) %>%
        rbindlist() %>%
        dcast(., featureID ~ sample, value.var = "count", fun.aggregate = mean, fill = 0) %>%
        select(featureID, all_of(sampleNames)) %>%
        tibble::column_to_rownames(var = "featureID")
}

generateColDataDF <- function(wholeFilesPath, sampleNames, condition, replicate, ...) {
    .check_argTypes(
        args = list(wholeFilesPath, sampleNames, condition, replicate),
        expected_types = c("character", "character", "character", "character")
    )
    .check_lengths(wholeFilesPath, sampleNames, condition, replicate)

    coldata <- data.frame(
        sampleLocation = wholeFilesPath,
        sampleName = sampleNames,
        condition = factor(condition, levels = unique(condition)),
        replicate = factor(replicate, levels = unique(replicate)),
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

generateComparisonDF <- function(pairList, join = TRUE) {
    if (!is.list(pairList) || !all(sapply(pairList, function(x) length(x) == 2 && all(sapply(x, is.character))))) {
        stop("Input must be a list of pairs: list(c('treat', 'control'))")
    }

    treat <- sapply(pairList, function(x) x[1])
    control <- sapply(pairList, function(x) x[2])

    comparisonDF <- data.frame(Treat = treat, Control = control, stringsAsFactors = FALSE)
    if (join) {
        comparisonDF$comparison <- paste0(treat, "_vs_", control)
    }
    comparisonDF
}




### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### utils
###
background <- (theme_bw() + theme(
    plot.title = element_text(size = 16),
    axis.title.x = element_text(size = 14, vjust = 0.5, hjust = 0.5, colour = "black"),
    axis.title.y = element_text(size = 14, vjust = 0.5, hjust = 0.5, colour = "black"),
    axis.text.x = element_text(size = 11, vjust = 0.5, hjust = 0.5, colour = "black"),
    axis.text.y = element_text(size = 11, vjust = 0.5, hjust = 0.5, colour = "black"),
    axis.line = element_line(linewidth = 0.5), legend.background = element_blank(),
    strip.text = element_text(size = 14, vjust = 0.5, hjust = 0.5, colour = "black"),
    strip.background = element_rect(fill = NA, colour = NA),
    panel.grid = element_blank()
))
