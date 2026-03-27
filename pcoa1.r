# -*- coding: utf-8 -*-
# 1. 安装并加载必要的包
# 如果你还没有安装这些包，请取消下面几行的注释并运行一次：
# install.packages(c("tidyverse", "vegan", "ggforce", "reshape2"))

library(tidyverse) # 数据处理
library(vegan)     # 多元统计 (vegdist, cmdscale)
library(ggforce)   # 用于画置信圈 (geom_mark_ellipse)
library(reshape2) # 数据重塑 (melt)

# 2. 数据读取与预处理
cat("正在读取数据...\n")

# 读取 Excel 文件 (请确保文件在当前工作目录下，或填写完整路径)
# 注意：你的文件名是 2pcoa(1).xlsx，R中文件名如果有特殊字符或括号，建议改名或用引号括起来
# 这里假设文件在当前目录下
raw_data <- readxl::read_excel("2pcoa(2).xlsx", sheet = 1)

# --- 数据清洗部分 ---
# 2.1 提取 Taxonomy 作为行名，并删除第一列
## 处理重复行名：自动加后缀保证唯一
raw_names <- raw_data[[1]]
if(any(duplicated(raw_names))) {
  raw_names <- make.unique(as.character(raw_names))
  cat("检测到重复行名，已自动修正为唯一。\n")
}
rownames(raw_data) <- raw_names
data_matrix <- raw_data %>% 
  select(-1) # 移除第一列 (Taxonomy)，行名已在前面设置，无需再 column_to_rownames

# 2.2 转置数据 (因为 PCoA 需要 样本 x 物种 矩阵，行是样本，列是物种/OTU)
# 注意：你的数据结构是 OTU x Sample，所以我们需要转置
if (!is.matrix(data_matrix)) {
  data_matrix <- as.matrix(data_matrix)
}
data_t <- t(data_matrix) # 转置后，行是样本，列是 OTU

# 3. 样本名称解析与分组 (关键步骤)
# 你的样本名如 S1Y1.1, M1Y1.1，我们提取 S1, M1 作为组别
cat("正在解析样本分组...\n")

# 提取样本名 (行名)
sample_names <- rownames(data_t)

# 使用正则表达式提取前缀 (匹配开头的非数字字符+第一个数字)
# 例如: S1Y1.1 -> S1; M1Y1.1 -> M1
group_matches <- str_match(sample_names, pattern = "^([A-Za-z]+\\d+)")

# 提取匹配的第一组 (如果有匹配)
groups <- ifelse(is.na(group_matches[, 2]), "Unknown", group_matches[, 2])

# 创建样本元数据 (Sample Metadata)
meta_data <- data.frame(
  Sample = sample_names,
  Group = factor(groups), # 转为因子，用于绘图颜色
  stringsAsFactors = FALSE
)

# 4. 计算距离矩阵 (Bray-Curtis)
cat("正在计算距离矩阵...\n")
# 移除可能存在的零和行
data_t <- data_t[rowSums(data_t) > 0, ]
# 计算距离
dist_matrix <- vegdist(data_t, method = "bray")

# 5. PCoA 分析 (经典多维尺度分析)
cat("正在执行 PCoA 分析...\n")
pcoa_result <- cmdscale(dist_matrix, k = 2, eig = TRUE) # k=2 表示提取前2个轴

# 提取特征值信息 (用于计算解释率)
eigenvalues <- pcoa_result$eig
total_var <- sum(eigenvalues[eigenvalues > 0])
pct <- round(eigenvalues[1:2] / total_var * 100, 2) # 前两个轴的百分比

# 提取坐标点数据
points <- pcoa_result$points
pcoa_df <- data.frame(
  PCoA1 = points[, 1],
  PCoA2 = points[, 2],
  Sample = rownames(points),
  stringsAsFactors = FALSE
)

# 6. 合并数据以便绘图
final_df <- merge(pcoa_df, meta_data, by.x = "Sample", by.y = "Sample")

# 7. 绘图 (使用 ggforce 画置信圈)
cat("正在绘制带置信圈的 PCoA 图...\n")

p <- ggplot(final_df, aes(x = PCoA1, y = PCoA2)) +
  
  # --- 核心：使用 ggforce 画置信圈 ---
  # geom_mark_ellipse 比 stat_ellipse 更适合小样本，且有膨胀效果
  geom_mark_ellipse(aes(fill = Group, color = Group, group = Group),
                    alpha = 0.1,           # 圈内透明度
                    expand = unit(1, "mm"), # 边缘向外膨胀 1mm，防止圈太紧
                    show.legend = FALSE) +  # 不显示圈的图例 (避免重复)
  
  # --- 散点图 ---
  geom_point(aes(color = Group), size = 4) +
  
  # --- 样本标签 (可选，样本太多可能会重叠，建议注释掉或调整) ---
  # geom_text(aes(label = Sample, color = Group), vjust = -0.5, size = 3) +
  
  # --- 坐标轴主题 ---
  scale_color_brewer(palette = "Set1") + # 点的颜色
  scale_fill_brewer(palette = "Set1") +  # 圈的颜色 (guide = "none" 已在上面隐藏)
  theme_minimal() +
  labs(
    title = "PCoA Plot with Confidence Ellipses (Bray-Curtis)",
    subtitle = "Grouped by Sample Prefix (e.g., S1, M1, L1)",
    x = paste0("PCoA 1 (", pct[1], "%)"),
    y = paste0("PCoA 2 (", pct[2], "%)"),
    color = "Group"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    axis.title = element_text(size = 12),
    legend.position = "right",
    # 加粗坐标轴线条
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_blank()
  ) +
  # 保持纵横比一致
  coord_fixed()

print(p)

# 8. 保存图片
ggsave("PCoA_With_Ellipse.png", plot = p, width = 12, height = 10, dpi = 300)
cat("分析完成！图片已保存为 \"D:\\R\\love-eat\\PCoA_With_Ellipse.png\"\n")