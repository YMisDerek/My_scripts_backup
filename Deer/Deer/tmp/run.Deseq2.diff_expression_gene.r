rm(list = ls())

### 传参
args <- commandArgs(trailingOnly = TRUE)
countsFileName <- args[1]               # matrix filename
groupsToExtract <- args[2]              # "JQ1/WT; " (from condition colunm)
filter <- as.numeric(args[3])
cutoff <- as.numeric(args[4])
refSpecies <- as.character(args[5])

### 参数处理
signif <- 0.05
coldataFileName <- sub("matrix", "coldata", countsFileName)
groups.AvsB <- strsplit(groupsToExtract, ";\\s*")[[1]]
# outputFileName <- sub("matrix", paste0("Deseq2_res",'.log2FC_',cutoff,'.filt_',filter), countsFileName)
outputFileName <- sub("counts.matrix", paste0("DEG",'.log2FC_',cutoff,'.filt_',filter), countsFileName)

### 检查、打印参数
if (file.exists(countsFileName) && file.exists(coldataFileName)) {
    message("========================DEseq2==========================")
    message("\n*****Working Dir")
    message(getwd())
    message("\n*****INPUT")
    message("counts: ", "\t", countsFileName)
    message("colData:", "\t", coldataFileName)
    message("\n*****ARGS")
    message("filter out counts less than: ", filter)
    message("significant cutoff:          ", signif)
    message("log2foldchange cutoff:       ", cutoff)
    message(length(groups.AvsB), " groups to compare:", "\t", paste0(groups.AvsB, "\t"))
    message("\n*****OUTPUT")
    message("result: ", "\t", paste0(outputFileName, "*"))
    message("========================================================\n")
} else {
    stop("counts file not exist;\ncheck working Directory.")
}



### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
# ***************************
#     差异分析
# ***************************
suppressWarnings(suppressMessages(library(DESeq2)))
suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(tidyr)))
suppressWarnings(suppressMessages(library(tibble)))
suppressWarnings(suppressMessages(library(purrr)))
suppressWarnings(suppressMessages(library(stringr)))
suppressWarnings(suppressMessages(library(ggplot2)))
suppressWarnings(suppressMessages(library(lay)))
suppressWarnings(suppressMessages(library(ggrastr)))
# source("~/.R/.self_priority.R")
load("~/.R/.self_theme.Rdata")

###STEP 1; DE analysis
colData <- read.table(coldataFileName, header = T, stringsAsFactors = F)
colData$condition <- factor(colData$condition, levels = unique(colData$condition))
colData$colNames <- factor(colData$colNames, levels = unique(colData$colNames))
rawCounts <- data.table::fread(countsFileName, data.table = F) %>%
    column_to_rownames(var = colnames(.)[1])
counts <- rawCounts[rowSums(rawCounts) >= filter, ]
colnames(counts) <- colData$colNames


dds <- DESeqDataSetFromMatrix(countData = counts, colData = colData, design = ~condition)
dds <- DESeq(dds, fitType = "local")
rld <- rlog(dds, blind = FALSE)
counts.rlog <- assay(rld) %>% as_tibble(., rownames = "ID")



###STEP 2; 整体统计，可视化
pdf(paste0(outputFileName, '.pdf'), width = 15)

fig.pca.all <-
    plotPCA(rld) +
    background +
    theme(legend.position = "none", aspect.ratio = 1) +
    coord_cartesian() +
    ggrepel::geom_label_repel(aes(label = name), size = 8 / .pt)
fig.pca <-
    plotPCA(rld, returnData = T) %>%
    filter(condition %in% unique(unlist(strsplit(groups.AvsB, '/')))) %>%
    ggplot(aes(PC1, PC2, color = condition)) +
    geom_point() +
    background +
    theme(legend.position = "none", aspect.ratio = 1) +
    coord_cartesian() +
    ggrepel::geom_label_repel(aes(label = name), size = 10 / .pt)
pca.variance <- summary(prcomp(assay(rld)))$importance["Proportion of Variance", c("PC1", "PC2")]
fig.pca <-
    fig.pca +
    labs(x = paste0('PC1: ', round(pca.variance['PC1']*100), '% variance'),
         y = paste0('PC2: ', round(pca.variance['PC2']*100), '% variance'))
fig.count <-
    pivot_longer(counts.rlog, -ID) %>%
    ggplot(aes(x = factor(name, levels = colData$colNames), y = value)) +
    geom_boxplot(notch = T, fill = "grey") +
    labs(x = NULL, y = "normalized signal") +
    background +
    theme(axis.text.x = element_text(angle = 90, hjust = 1), aspect.ratio = 1)
samplesCorrDF <- 
    Hmisc::rcorr(as.matrix(counts.rlog[, -1]))$r %>% 
    as.data.frame() %>% rownames_to_column() %>%
    pivot_longer(cols = -1) %>% 
    rename_with(~c("name1", "name2", "rho"))
fig.corr <- 
    ggplot(samplesCorrDF) +
    geom_tile(aes(factor(name1, levels = colData$colNames), 
                  factor(name2, levels = colData$colNames), 
                  fill = rho)) +
    geom_text(aes(factor(name1, levels = colData$colNames), 
                  factor(name2, levels = colData$colNames), 
                  label = round(rho, 2)), size = 2) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1),
          axis.text = element_text(color = "black"),
          axis.title = element_blank(),
          aspect.ratio = 1) +
    scale_fill_gradient(low = "white", high = "#1b98e0")

print(ggpubr::ggarrange(fig.pca.all, fig.pca, nrow = 1, ncol = 2, align = "h"))
print(ggpubr::ggarrange(fig.count, fig.corr, nrow = 1, ncol = 2, align = "h"))



###STEP 3; 提取组间比较结果，并统计绘图
# 准备提取的组合信息
groups.AvsB.all.df <-
    expand.grid(colData$condition, colData$condition) %>%
    subset(Var1 != Var2) %>%
    distinct() %>%
    mutate(AvsB = paste0(Var1,"/",Var2))
if (all(groups.AvsB %in% groups.AvsB.all.df$AvsB)) {
    groups.AvsB.all.df <-
        groups.AvsB.all.df %>%
        arrange(match(AvsB, groups.AvsB)) %>%
        filter(AvsB %in% groups.AvsB)
    message('\n\n\n*****GROUP TO VS')
    print(groups.AvsB.all.df)
} else {
    message("# valid groups\n\n", paste(groups.AvsB.all.df$AvsB, "\n"))
    stop("<Input group AvsB> NOT in <coldata conditions>")
}

# 每组统计绘制火山图、散点图
for(i in 1:nrow(groups.AvsB.all.df)){
    conA <- groups.AvsB.all.df$Var1[i] %>% as.character()
    conB <- groups.AvsB.all.df$Var2[i] %>% as.character()
    message(conA, ' vs ', conB)

    # 提取组结果
    deseq2.res.AvsB <- 
        results(dds, contrast = c("condition", conA, conB)) %>%
        as_tibble(rownames = "ID") %>%
        mutate(type = case_when(padj < signif & log2FoldChange > cutoff ~ "up",
                                padj < signif & log2FoldChange < cutoff * (-1) ~ "down",
                                TRUE ~ "None")) %>%
        dplyr::select(1, 3, 7, 8) %>%
        rename_with(~ c("ID", "log2FC", "padj", "deType"))

    # 统计绘图
    changeNum <-
        dplyr::count(deseq2.res.AvsB, deType, name = "Freq") %>%
        mutate(tag = paste0(deType, ': ', Freq))
    xlab_range <- deseq2.res.AvsB$log2FC %>% abs() %>% max() %>% ceiling()
    ylab_range <- deseq2.res.AvsB$padj %>% -log10(.) %>% max() %>% ceiling()

    volc <-
        ggplot(deseq2.res.AvsB) +
        geom_point_rast( aes(x = log2FC, y = -log10(padj), color = deType), size = 0.3) +
        geom_vline(xintercept = c(0), color = "#990000", linetype = "dashed") +
        scale_color_manual(values = c('down'= "#4d97cd", 'None'="#dee2e6", 'up'="#c74546"),
                           labels = setNames(changeNum$tag, changeNum$deType),
                           name = paste0(conA, '\nvs\n', conB)) +
        scale_x_continuous(limits = c(-xlab_range, xlab_range),
                           breaks = seq(-xlab_range, xlab_range, by = 1)) +
        background +
        theme(aspect.ratio = 1) +
        guides(color = guide_legend(override.aes = list(size = 5)))

    sca <-
        dplyr::mutate(counts.rlog, 
                      treat_norm = lay(across(contains(conA)), mean), 
                      control_norm = lay(across(contains(conB)), mean)) %>%
        dplyr::select(ID, treat_norm, control_norm) %>%
        dplyr::left_join(., dplyr::select(deseq2.res.AvsB, ID, deType), by = "ID") %>%
        ggplot(aes(x = control_norm, y = treat_norm, color = deType)) +
        geom_point_rast(size = 0.2, alpha = 0.7) +
        geom_abline(slope = 1, intercept = 0, lty = 2) +
        labs(x = paste0("Intensity of ", conB), y = paste0("Intensity of ", conA), color = "deType") +
        scale_color_manual(values = c('down'= "#4d97cd", 'None'="#dee2e6", 'up'="#c74546")) +
        scale_x_continuous(breaks = c(0, 3, 6, 9, 12, 15), limits = c(0, 15)) +
        scale_y_continuous(breaks = c(0, 3, 6, 9, 12, 15), limits = c(0, 15)) +
        background +
        guides(color = guide_legend(override.aes = list(size = 5))) +
        theme(legend.position = "none", aspect.ratio = 1)
    
    print(ggpubr::ggarrange(volc, sca, nrow = 1, ncol = 2, align = 'h', common.legend = T, legend = "right"))

    # 保存结果
    names(deseq2.res.AvsB) <- c("ID", paste(c("log2FC", "padj", "deType"), paste0(conA,'_vs_', conB), sep = "."))
    if(exists("deseq2.res")){
        deseq2.res <- left_join(deseq2.res, deseq2.res.AvsB, by = 'ID')
    }else{
        deseq2.res <- deseq2.res.AvsB
    }
}
dev.off()

### save res
names(counts.rlog) <- paste0(names(counts.rlog), '.norm')
deseq2.res <- cbind(deseq2.res, counts, counts.rlog[, -1])

if(refSpecies == "human"){
    deseq2.res <- 
        dplyr::left_join(deseq2.res, annotables::grch38, by = c("ID" = "ensgene")) %>% 
        dplyr::group_by(ID) %>% 
        dplyr::slice(1) %>% 
        dplyr::ungroup()
}else if(refSpecies == "mouse"){
    deseq2.res <- 
        dplyr::left_join(deseq2.res, annotables::grcm38, by = c("ID" = "ensgene")) %>% 
        dplyr::group_by(ID) %>% 
        dplyr::slice(1) %>% 
        dplyr::ungroup()
}


write.table(deseq2.res, paste0(outputFileName, ".xls"), quote = F, sep = "\t", row.names = F, col.names = T)
resList= list(dds=dds, colData=colData, groups.AvsB=groups.AvsB, deseq2.res=deseq2.res)
saveRDS(resList, file = paste0(outputFileName, '.rds'))
