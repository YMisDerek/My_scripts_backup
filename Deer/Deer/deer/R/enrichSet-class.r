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