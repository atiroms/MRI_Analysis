library("ggpubr")


ggscatter(e, x = "dim1", y = "age", add = "reg.line", conf.int = TRUE,cor.coef = TRUE, cor.method = "pearson")