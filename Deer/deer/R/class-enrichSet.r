### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### class
###

setClass("enrichSet",
    slots = list(
        enrichResults = "list",
        enrich.res = "data.frame",
        enrich.image = "list"
    )
)



### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### show
###

setMethod("show", "enrichSet", function(object) {
    cat("\n..@enrichSet\n")
    cat("class:", class(object), "\n")
    # cat("Diff result: ", dim(object@DE.result), "\n")

    # coolcat("Diff statistic(%d): %s\n", names(object@DE.stat))
    # coolcat("Diff image(%d): %s\n", names(object@DE.image))
})



### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### method
###

# combine all enrichment to DF for PyComplexHeatmap plot
setGeneric("exportDF", signature = "obj", function(obj) standardGeneric("exportDF"))

.export_DF_PyComplexHeatmap = function(obj) {

}

setMethod("exportDF", "enrichSet", .export_DF_PyComplexHeatmap)