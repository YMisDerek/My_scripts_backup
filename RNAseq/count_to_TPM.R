# ************************
# count to tpm
# ************************
gtf.mm10 <- '/data/public/Gene.annotation/mm10/Mus_musculus.GRCm38.98.gtf'
# txdb.mm10 <- makeTxDbFromGFF(gtf.mm10)
# saveDb(txdb.mm10, 'mm10.grcm38.98.gtf.txdb.qslite')
txdb.mm10 <- loadDb('mm10.grcm38.98.gtf.txdb.qslite')

normalizeGeneCounts <- function(counts, TxDb, method) {
    require("GenomicFeatures")
    
    # len <- width(genes(TxDb))/1e3
    len <- sum(width(exonsBy(TxDb, by = 'gene')))/1e3
    
    if (method == "CPM") {
        normalized <-
            apply(counts,2, function(x) (x/sum(x)) * 1e6)
    } else if (method == "RPKM") {
        normalized <-
            t(apply(counts, 1, "/", colSums(counts) / 1e6)) / len
    } else if (method == "TPM") {
        x <- counts / len
        normalized <- t( t(x)*1e6 / colSums(x) )
    } else {
        stop("method must be CPM, RPKM or TPM")
    }
    return(as.data.frame(normalized))
}
