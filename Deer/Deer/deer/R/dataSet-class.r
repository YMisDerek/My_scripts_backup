### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### dataSet object
###
### 针对不同数据类型创建dataSet的功能还可以改进
### 暂时依照assay.type类型判断，且最终均转为RangedSummarizedExperiment

setClass("dataSet",
    contains = "SummarizedExperiment",
    slots = list(
        project = "character_OR_NULL",
        organism = "character_OR_NULL"
    )
)

createDataSet <- function(countMatrix, ColDataDF, RowDataDF, comparisonDF,
                          assay.type, assay.design, min.count, DE.foldchange, DE.qvalue) {
    if (assay.type == "RNA") {
        .check_argTypes(list(RowDataDF), c("data.frame"))

        dataSet <- SummarizedExperiment(
            assays = list(counts = countMatrix),
            colData = ColDataDF,
            rowData = RowDataDF,
            metadata = list(
                "assay.type" = assay.type,
                "assay.design" = ColDataDF,
                "comparisonDF" = comparisonDF,
                "min.count" = min.count,
                "DE.foldchange" = DE.foldchange,
                "DE.qvalue" = DE.qvalue,
                "active.features" = rownames(countMatrix)[rowSums(countMatrix) > min.count],
                "active.samples" = ColDataDF$sampleName
            )
        )
        dataSet <- as(dataSet, "RangedSummarizedExperiment")
    } else if (assay.type == 'ChIP') {
        .check_argTypes(list(RowDataDF), c("GRangesList"))
        # ...
    }
    return(dataSet)
}
