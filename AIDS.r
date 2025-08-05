x <- read.table("x.txt", header = TRUE)
predelta <- read.table("predelta.txt", header = TRUE)
colnames(predelta) <- "CATE"
df <- cbind(x, CATE = predelta$CATE)

# 加载绘图包
library(ggplot2)
library(patchwork)

# 图 1
p1 <- ggplot(df, aes(x = Age, y = CATE, color = as.factor(Tt))) +
    geom_point(alpha = 0.7) +
    labs(title = "CATE vs Age", x = "Age", y = "Estimated CATE", color = "Treatment") +
    theme_minimal()

# 图 2
p2 <- ggplot(df, aes(x = CD4_0, y = CATE, color = as.factor(Tt))) +
    geom_point(alpha = 0.7) +
    labs(title = "CATE vs Baseline CD4 Count", x = "CD4_0", y = "Estimated CATE", color = "Treatment") +
    theme_minimal()

# 合并图 + 图例置底
(p1 / p2) + 
    plot_layout(guides = "collect") & 
    theme(legend.position = "bottom")
