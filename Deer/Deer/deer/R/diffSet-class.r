### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### class
###

setClass("diffSet",
    slots = list(
        dds = "DESeqDataSet",
        DE.result = "data.frame",
        DE.stat = "list",
        DE.image = "list"
    )
)




### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### show
###

setMethod(
    "show", "diffSet",
    function(object) {
        cat("class:", class(object), "\n")
        cat("Diff result: ", dim(object@DE.result), "\n")
        cat("Diff statistic(): ", names(object@DE.stat), "\n")
        cat("Diff image(): ", names(object@DE.image), "\n")

        # cat("\n...@dds...\n")
        # print(object@dds)

        ## DE.result
    }
)
