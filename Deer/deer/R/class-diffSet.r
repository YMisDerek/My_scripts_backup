### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### class
###

# setClass("diffSet",
#     contains = "DESeqDataSet",
#     slots = list(
#         DE.result = "data.frame",
#         DE.stat = "list",
#         DE.image = "list"
#     )
# )
setClass("diffSet",
    slots = list(
        DE.dds = "ANY",
        DE.stat = "data.frame",
        DE.image = "list"
    )
)



### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### show
###

setMethod("show", "diffSet", function(object) {
    cat("\n...@diffSet\n")

    cat("......Diff stat: \n")
    if (nrow(object@DE.stat) > 0) print(object@DE.stat)
    coolcat("......Diff image(%d): %s\n", names(object@DE.image))
})
