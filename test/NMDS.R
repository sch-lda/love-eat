rm(list = ls())
library(openxlsx)
library(ggplot2)
library(ggpubr)
library(ggforce)
library(cowplot)
library(patchwork)
library(vegan)

## 1) 读入数据 ---------------------------------------------------------------
otu_data <- read.xlsx("./Mun-number.xlsx",
                      sheet = 1, rowNames = TRUE, colNames = TRUE)
group_data <- read.xlsx("./Reads-number.xlsx",
                        sheet = 2)

## 2) 距离矩阵（Bray-Curtis）
dist_matrix <- vegdist(t(otu_data), method = "bray")

## 3) NMDS 拟合 --------------------------------------------------------------
set.seed(123)  # 可复现
nmds_result <- metaMDS(dist_matrix,
                       k = 2,
                       trymax = 200,         # 提高迭代上限更稳妥
                       autotransform = FALSE,# 我们已提供距离矩阵
                       wascores = TRUE,
                       trace = FALSE)

## 4) 提取坐标与合并分组信息 -------------------------------------------------
# metaMDS 的站点坐标（样本）
nmds_scores <- as.data.frame(scores(nmds_result, display = "sites"))
nmds_scores$Site <- rownames(nmds_scores)

# 与分组信息按样本名安全合并（假定 group_data 有列“Site”与“Group”等）
plot_data <- merge(nmds_scores, group_data, by = "Site", all.x = TRUE)

## 5) 关键统计量：stress、R²、逐样本 goodness -----------------------------------
# stress（Kruskal）
nmds_stress <- nmds_result$stress

# 观测差异 vs 排序距离的相关决定系数（常用总结指标）
obs_d <- as.vector(dist_matrix)
ord_d <- as.vector(dist(nmds_scores[, c("NMDS1", "NMDS2")]))
nmds_r2 <- cor(obs_d, ord_d, method = "pearson")^2

# 逐样本 goodness（值越小拟合越好）
gof <- goodness(nmds_result)              # 命名向量
gof_df <- data.frame(Site = names(gof), goodness = as.numeric(gof))

# 输出到控制台
cat("NMDS stress =", round(nmds_stress, 4), "\n")
cat("NMDS fit R^2 (dist vs ord-dist) =", round(nmds_r2, 4), "\n")
head(gof_df)

# 如需导出逐样本 goodness：
# write.csv(gof_df, "NMDS_goodness_by_sample.csv", row.names = FALSE)

## 6) PERMANOVA（与 PCoA 一致的检验方式） ------------------------------------
# 注意：adonis2 使用“原始矩阵”的样本元数据（group_data），公式里因子名改成你的列名
dunediv <- adonis2(t(otu_data) ~ Group, data = group_data,
                   permutations = 999, method = "bray")
print(dunediv)

## 7) 可视化（ggplot + geom_mark_ellipse） ------------------------------------
# 图标题里标注 stress；坐标轴名为 NMDS1/NMDS2
subtitle_txt <- paste0("Stress = ", round(nmds_stress, 3),
                       " | Fit R^2 = ", round(nmds_r2, 3))

p_nmds <- ggplot() +
      geom_point(data = plot_data, size = 3,
                 aes(x = NMDS1, y = NMDS2, fill = Group, colour = Group)) +
      # 相切样式的分组包络：把 expand 设小一点更靠近点
      geom_mark_ellipse(data = plot_data,
                        aes(x = NMDS1, y = NMDS2, fill = Group, colour = Group),
                        expand = unit(1.5, "mm"),
                        show.legend = FALSE) +
      geom_vline(xintercept = 0, colour = "black", size = 0.5, linetype = 2) +
      geom_hline(yintercept = 0, colour = "black", size = 0.5, linetype = 2) +
      theme_bw() +
      labs(title = "NMDS (Bray-Curtis)",
           subtitle = subtitle_txt,
           x = "NMDS1", y = "NMDS2", fill = "", colour = "") +
      theme(text = element_text(family = "serif"),
            axis.title = element_text(size = 18),
            axis.text = element_text(size = 16),
            plot.title = element_text(size = 20),
            legend.title = element_blank(),
            panel.grid = element_blank(),
            legend.key.size = unit(10, "pt"),
            legend.background = element_blank()) +
      guides(fill = "none", colour = "none")

print(p_nmds)

## 8) （可选）Shepard 图检视非度量拟合 ---------------------------------------
# 会画出观测差异 vs 排序距离的散点与光滑曲线，便于评估非线性与离群点
# stressplot(nmds_result)
