####qqdeg
#自定义函数#
#' Title
#'
#' @param file file in document
#' @param object_type gene data or te data
#' @param group1 treat group
#' @param group2 control group
#' @param fc_threshold foldchange value
#'
#' @return deg results and enrich results
#' @export
#'
#' @examples
#' result <- qqdeg("rlim.xlsx","gene","male-ko","male-wt",fc_threshold = 1.5)
qqdeg <- function(file, object_type, group1, group2, fc_threshold = 1.5) {
  # 加载必要的包，并抑制启动消息
  suppressPackageStartupMessages({
    library(DESeq2)
    library(ggplot2)
    library(ggpubr)
    library(openxlsx)
    library(clusterProfiler)
    library(org.Mm.eg.db)
    library(ggthemes)
    library(ggrepel)
    library(gmodels)
    library(matrixStats)
    library(DOSE)
    library(topGO)
    library(pathview)
  })

  # 创建输出文件夹
  output_dir <- paste0(group1, "_vs_", group2, "_", object_type)
  dir.create(output_dir, showWarnings = FALSE)

  # 读取数据
  data <- read.xlsx(file, colNames = TRUE, rowNames = FALSE)

  # 根据对象类型提取数据
  if (object_type == "gene") {
    obj_data <- data[grep(data[,1], pattern = ":", invert = TRUE), ]
  } else if (object_type == "te") {
    obj_data <- data[grep(data[,1], pattern = ":", invert = FALSE), ]
  } else {
    stop("Invalid object_type. Must be 'gene' or 'te'.")
  }

  row.names(obj_data) <- obj_data$Gene
  obj_data <- obj_data[, -1]

  # 准备分组信息
  coldata <- data.frame(condition = gsub("-\\d+$", "", colnames(obj_data)), row.names = colnames(obj_data))

  # 创建DESeq2数据集
  dds <- DESeqDataSetFromMatrix(countData = obj_data, colData = coldata, design = ~condition)
  keep <- rowSums(counts(dds) >= 10) >= 2
  dds <- dds[keep, ]

  # 标准化和差异表达分析
  dds <- DESeq(dds)
  res <- results(dds, contrast = c("condition", group1, group2))

  # 保存差异分析结果
  resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized = TRUE)), by = "row.names", sort = FALSE)
  write.xlsx(resdata, file = file.path(output_dir, paste0(object_type, "_", group1, "_vs_", group2, ".xlsx")))

  # 提取显著差异基因
  diff_genes <- subset(res, padj < 0.05 & (log2FoldChange > log2(fc_threshold) | log2FoldChange < -log2(fc_threshold)))
  write.csv(diff_genes, file = file.path(output_dir, paste0("diff_", object_type, "_", group1, "_vs_", group2, ".csv")))

  # 提取上调和下调基因
  up_genes <- subset(res, padj < 0.05 & log2FoldChange > log2(fc_threshold))
  down_genes <- subset(res, padj < 0.05 & log2FoldChange < -log2(fc_threshold))
  write.csv(up_genes, file = file.path(output_dir, paste0(group1, "_vs_", group2, "_up-", object_type, ".csv")))
  write.csv(down_genes, file = file.path(output_dir, paste0(group1, "_vs_", group2, "_down-", object_type, ".csv")))

  # 绘制火山图
  resdata$logP <- -log10(resdata$padj)
  resdata$GROUP <- "not significant"
  resdata$GROUP[(resdata$padj < 0.05) & (resdata$log2FoldChange > log2(fc_threshold))] <- "up-regulated"
  resdata$GROUP[(resdata$padj < 0.05) & (resdata$log2FoldChange < -log2(fc_threshold))] <- "down-regulated"
  resdata$Lable <- ""

  # 去除缺失值
  resdata <- na.omit(resdata)
  resdata <- resdata[order(resdata$padj), ]
  top_up_genes <- head(resdata$Row.names[resdata$GROUP == "up-regulated"], 10)
  top_down_genes <- head(resdata$Row.names[resdata$GROUP == "down-regulated"], 10)
  top_genes <- c(as.character(top_up_genes), as.character(top_down_genes))
  resdata$Lable[match(top_genes, resdata$Row.names)] <- top_genes

  palette <- c("#2f5688", "#BBBBBB", "#CC0000")
  group_colors <- c("down-regulated" = palette[1], "not significant" = palette[2], "up-regulated" = palette[3])

  volcano_plot <- ggplot(resdata, aes(x = log2FoldChange, y = logP)) +
    geom_point(aes(color = GROUP), size = 1) +
    geom_text_repel(aes(label = Lable), size = 3, box.padding = 0.5, max.overlaps = Inf) +
    scale_color_manual(values = group_colors) +
    theme_classic() +
    geom_hline(yintercept = 1.30, linetype = "dashed") +
    geom_vline(xintercept = c(-log2(fc_threshold), log2(fc_threshold)), linetype = "dashed") +
    expand_limits(y = 0) +
    xlab(bquote(log[2]*FoldChange * "(" ~ .(group1) ~ "/" ~ .(group2) ~ ")")) +
    ylab(expression(-log[10](Adjusted ~ P-value)))

  # 显示火山图
  print(volcano_plot)

  # 保存火山图
  ggsave(filename = file.path(output_dir, paste0("Volcano_plot_", object_type, "_", group1, "_vs_", group2, ".png")), plot = volcano_plot, limitsize = FALSE)

  # GO富集分析
  go <- NULL
  kk <- NULL
  if (object_type == "gene") {
    up <- read.csv(file.path(output_dir, paste0(group1, "_vs_", group2, "_up-", object_type, ".csv")))
    down <- read.csv(file.path(output_dir, paste0(group1, "_vs_", group2, "_down-", object_type, ".csv")))

    up_genes <- bitr(up$X, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db, drop = FALSE)
    down_genes <- bitr(down$X, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db, drop = FALSE)

    if (nrow(up_genes) == 0 | nrow(down_genes) == 0) {
      cat("Warning: No genes mapped to ENTREZ IDs. Skipping GO and KEGG analysis.\n")
    } else {
      up.go_all <- enrichGO(gene = up_genes$ENTREZID, OrgDb = org.Mm.eg.db, keyType = 'ENTREZID', ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1)
      down.go_all <- enrichGO(gene = down_genes$ENTREZID, OrgDb = org.Mm.eg.db, keyType = 'ENTREZID', ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1)

      # 处理GO结果
      go.up <- up.go_all@result
      go.up <- go.up[order(go.up$p.adjust), ]
      go.up$value <- -log10(go.up$p.adjust)

      go.down <- down.go_all@result
      go.down <- go.down[order(go.down$p.adjust), ]
      go.down$new <- -log10(go.down$p.adjust)
      go.down$value <- -(go.down$new)
      go.down <- go.down[ , -which(colnames(go.down) %in% "new")]

      # 合并显著的上下调
      go.up <- head(go.up, 10)
      go.down <- head(go.down, 10)
      go <- rbind(go.up, go.down)
      go <- go[go$value > 1.30103 | go$value < -1.30103, ]

      ## 调整因子水平
      go <- go[order(go$value),]
      # 计算每个 Description 的频率
      desc_counts <- table(go$Description)
      # 筛选出只出现一次的 Description
      unique_desc <- names(desc_counts[desc_counts == 1])
      # 根据唯一的 Description 过滤 go 数据框
      go <- go[go$Description %in% unique_desc, ]
      name <- go$Description
      go$Description <- factor(go$Description, levels = name)

      # 绘制GO富集图
      go_plot <- ggplot(go, aes(x = value, y = Description)) +
        xlab('Enrich value') +
        ylab("") +
        geom_bar(stat = "identity", aes(fill = ifelse(value < 0, "down", "up"))) +
        scale_fill_manual(name = "value", values = c("down" = "dodgerblue3", "up" = "firebrick3")) +
        theme_few() +
        ggtitle(paste(group1, "vs", group2, "GO BP Enrichment Analysis"))

      print(go_plot)
      ggsave(filename = file.path(output_dir, paste0("GO_BP_Enrichment_", object_type, "_", group1, "_vs_", group2, ".png")), plot = go_plot, limitsize = FALSE)

      #######kegg
      #上调kegg
      up.kk <- enrichKEGG(gene = up_genes$ENTREZID,organism = 'mmu',  pvalueCutoff = 1)
      #下调kegg
      down.kk <- enrichKEGG(gene = down_genes$ENTREZID,organism = 'mmu',  pvalueCutoff = 1)
      #画图
      upkegg <- up.kk@result
      up.kegg <- upkegg[order(upkegg$p.adjust),]
      up.kegg$value <- -log10(up.kegg$p.adjust)
      downkegg <- down.kk@result
      down.kegg <- downkegg[order(downkegg$p.adjust),]
      down.kegg$value <- log10(downkegg$p.adjust)
      ####合并显著的上下调
      kegg.up <-  head(up.kegg,10)
      kegg.down <- head(down.kegg,10)
      kegg <- rbind(kegg.up,kegg.down)
      kegg <- kegg[kegg$value > 1.30103 | kegg$value < -1.30103, ]

      kegg$Description <- gsub("-.*", "", kegg$Description)
      ####BAR
      ## 调整因子水平
      kegg <- kegg[order(kegg$value),]
      # 计算每个 Description 的频率
      desc_counts <- table(kegg$Description)
      # 筛选出只出现一次的 Description
      unique_desc <- names(desc_counts[desc_counts == 1])
      # 根据唯一的 Description 过滤 kegg 数据框
      kegg <- kegg[kegg$Description %in% unique_desc, ]
      name <- kegg$Description
      kegg$Description <- factor(kegg$Description, levels = name)
      #绘制富集图
      kegg_plot <- ggplot(kegg, aes(x = value, y = Description),col = col) +
        xlab('Enrich value') +
        ylab("") +
        geom_bar(stat = "identity", aes(fill = ifelse(value < 0, "down", "up"))) +
        scale_fill_manual(name = "value", values = c("down" ="dodgerblue3", "up" = "firebrick3")) +
        theme_few() +
        ggtitle(paste(group1, "vs", group2, "KEGG Enrichment Analysis"))
      print(kegg_plot)

      ggsave(filename = file.path(output_dir, paste0("KEGG_Enrichment_", object_type, "_", group1, "_vs_", group2, ".png")), plot = kegg_plot, limitsize = FALSE)
    }
  }

  # 将所有结果保存到列表
  result <- list(
    volcano_plot = volcano_plot,
    go_plot = if (exists("go_plot")) go_plot else NULL,
    kegg_plot = if (exists("kk_plot")) kegg_plot else NULL,
    diff_genes = file.path(output_dir, paste0("diff_", object_type, "_", group1, "_vs_", group2, ".csv")),
    up_genes = file.path(output_dir, paste0(group1, "_vs_", group2, "_up-", object_type, ".csv")),
    down_genes = file.path(output_dir, paste0(group1, "_vs_", group2, "_down-", object_type, ".csv")),
    obj_data = obj_data,
    coldata = coldata,
    res = res,
    resdata = resdata,
    go = if (exists("go")) go else NULL,
    kegg = if (exists("kegg")) kegg else NULL
  )
  # 打印完毕消息
  # 创建虚线
  line <- strrep("-", 50)

  # 输出虚线和文本
  cat(line, "\n")
  cat("数据分析完毕，已保存！\n")
  cat(line, "\n")
  return(result)

}
#
