EnrichAnalyze <- function(object, GOKEGG_Genes, GOKEGG_Keytype, GOKEGG_OrgDb, KEGG_Org, GSEA_species, GSEA_Genes, ...) {
    ### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ### Args example
    ###
    # deg.group <- split(deg$entrez, deg$deType.SWA_vs_NC)
    # input_GSEA_Genes.df <-
    #     filter(deg, padj.SWA_vs_NC != "NA", symbol != "NA") %>%
    #     arrange(desc(log2FC.SWA_vs_NC)) %>%
    #     distinct(symbol, .keep_all = T)
    # input_GSEA_Genes <- input_GSEA_Genes.df$log2FC.SWA_vs_NC
    # names(input_GSEA_Genes) <- input_GSEA_Genes.df$symbol
    # input_GOKEGG_Genes <- list("DEG.down" = unique(na.omit(deg.group$down)), "DEG.up" = unique(na.omit(deg.group$up)))
    # input_GOKEGG_Keytype <- "ENTREZID"
    # input_GOKEGG_OrgDb <- org.Mm.eg.db::org.Mm.eg.db
    # input_KEGG_Org <- "mmu"
    # input_GSEA_species <- "Mus musculus"

    # GO
    input_GOKEGG_Genes <- GOKEGG_Genes
    input_GOKEGG_Keytype <- GOKEGG_Keytype
    input_GOKEGG_OrgDb <- GOKEGG_OrgDb
    input_KEGG_Org <- KEGG_Org
    # GSEA
    input_GSEA_Genes <- GSEA_Genes
    input_GSEA_species <- GSEA_species

    enrich.go.res <- map_df(names(input_GOKEGG_Genes), function(name) {
        enrich.go <- enrichGO(OrgDb = input_GOKEGG_OrgDb, gene = input_GOKEGG_Genes[[name]], keyType = input_GOKEGG_Keytype, ont = "ALL", readable = T)
        enrich.go@result$DEG_type <- as.character(name)
        as.data.frame(enrich.go)
    })

    enrich.kegg.res <- map_df(names(input_GOKEGG_Genes), function(name) {
        enrich.kegg <- enrichKEGG(organism = input_KEGG_Org, gene = input_GOKEGG_Genes[[name]]) %>% setReadable(x = ., OrgDb = input_GOKEGG_OrgDb, keyType = "ENTREZID")
        enrich.kegg@result$DEG_type <- as.character(name)
        as.data.frame(enrich.kegg)
    })

    set.seed(2024)
    hs_hallmark_sets <- msigdbr::msigdbr(species = input_GSEA_species, category = "H")
    enrich.gsea <- clusterProfiler::GSEA(geneList = input_GSEA_Genes, seed = TRUE, TERM2GENE = dplyr::select(hs_hallmark_sets, gs_name, gene_symbol))
    enrich.gsea.df <- as.data.frame(enrich.gsea)

    object@enrichSet@enrichResults <- list("GO" = enrich.go.res, "KEGG" = enrich.kegg.res, "GSEA" = enrich.gsea.df)
}
