# 1. 清空环境并加载必要的包
rm(list = ls())
# 自动安装并加载必要的包
required_packages <- c("readxl", "vegan", "ggplot2", "dplyr")
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages, repos = "https://cloud.r-project.org")

library(readxl)   # 读取Excel
library(vegan)    # 计算距离和PCoA
library(ggplot2)  # 绘图
library(dplyr)    # 数据处理
library(ggforce)

# 2. 设置文件路径 (请将文件放在工作目录下，或写完整路径)
meta_file <- "./h1.xlsx"   # 包含 Site 和 Group 信息
otu_file  <- "./h2.xlsx"   # 包含 Taxonomy 和 丰度信息

# 3. 读取数据
cat("正在读取数据...\n")
meta_data <- read_excel(meta_file, sheet = 1)
otu_data  <- read_excel(otu_file,  sheet = 1)

# 4. 数据预处理 (关键步骤)
# 4.1 处理OTU表 (侯2.xlsx)
# 假设第一列是分类学信息，其余列是样本数据
otu_matrix <- as.matrix(otu_data[, -1]) # 去掉第一列Taxonomy
rownames(otu_matrix) <- otu_data[[1]]   # 将第一列作为行名

# 4.2 处理元数据 (侯1.xlsx)
# 提取样本ID (去除后缀，例如 S_1 -> S) 或者直接使用 Site
# 这里我们假设 侯1.xlsx 中的 "Site" 列 (S_1, S_2...) 对应 OTU表中的列名 (S_1, S_2...)
sample_ids <- meta_data$Site

# 4.3 校验：确保两份文件的样本顺序一致
# 提取OTU表中的列名
otu_colnames <- colnames(otu_matrix)

# 检查元数据中的Site是否都在OTU表中
missing <- setdiff(sample_ids, otu_colnames)
if(length(missing) > 0) {
  warning(paste("以下样本在OTU表中未找到:", paste(missing, collapse = ", ")))
}

# 重新排列OTU矩阵，使其列顺序与 meta_data 的行顺序一致
otu_sorted <- otu_matrix[, sample_ids, drop = FALSE]

# 5. 计算距离矩阵 (Bray-Curtis)
cat("正在计算Bray-Curtis距离...\n")
# 注意：数据应该是 变量(物种) x 样本 的矩阵，所以需要转置
dist_matrix <- vegdist(t(otu_sorted), method = "bray")

# 6. 执行 PCoA (主坐标分析)
cat("正在执行PCoA分析...\n")
pcoa_result <- capscale(dist_matrix ~ 1) # capscale 也可以做PCoA
# 或者使用 cmdscale: pcoa_result <- cmdscale(dist_matrix, k=10, eig=TRUE)

# 提取解释率 (用于坐标轴标签)
eig <- pcoa_result$CA$eig
pct <- round(eig / sum(eig) * 100, 2) # 前两个轴的百分比

# 提取样本坐标
coords <- as.data.frame(pcoa_result$CA$u)
# 明确命名坐标轴，以防 vegan 版本差异导致列名不匹配
colnames(coords)[1:2] <- c("PCoA1", "PCoA2")
coords$Sample <- rownames(coords)
# 反转PCoA1坐标轴（如需反转Y轴可同时加上下一行）
coords$PCoA2 <- -coords$PCoA2

# 7. 合并分组信息 (将 侯1.xlsx 的 Group 信息加入坐标表)
# 注意：coords 的行名是距离矩阵的名称，即 sample_ids
final_df <- coords %>%
  mutate(Sample = rownames(coords)) %>%
  left_join(meta_data %>% select(Sample = Site, Group), by = "Sample")

# 8. 绘图
cat("正在绘制PCoA图...\n")
p <- ggplot(final_df, aes(x = PCoA1, y = PCoA2, color = Group, label = Sample)) +
# --- 核心：使用 ggforce 画置信圈 ---
  # geom_mark_ellipse 比 stat_ellipse 更适合小样本，且有膨胀效果
  geom_mark_ellipse(aes(fill = Group, color = Group, group = Group),
                    alpha = 0.1,           # 圈内透明度
                    label.margin = margin(200, 200, 200, 200, "cm"), 
                    expand = unit(1, "mm"), # 边缘向外膨胀 1mm，防止圈太紧
                    show.legend = FALSE) +  # 不显示圈的图例 (避免重复)
  geom_point(size = 4, alpha = 0.8) +

  #geom_text(aes(label = Sample), vjust = -0.5, size = 3) + # 样本标签在点上方
  scale_color_viridis_d(option = "Set3") + # 使用更美观的颜色
  theme_minimal() +
  labs(
    title = "PCoA Plot (Bray-Curtis Distance)",
    x = paste0("PCoA 1 (", pct[1], "%)"),
    y = paste0("PCoA 2 (", pct[2], "%)"),
    color = "Group"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.title = element_text(size = 12),
    legend.position = "right"
  ) +
  # 添加网格线辅助观察
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray80") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray80") +
  expand_limits(x = range(final_df$PCoA1) + c(-0.2, 0.2),
                y = range(final_df$PCoA2) + c(-0.2, 0.2))

print(p)

# 9. (可选) 保存图片
ggsave("PCoA_Result.png", plot = p, width = 10, height = 8, dpi = 300)
cat("分析完成！图片已保存为 PCoA_Result.png\n")