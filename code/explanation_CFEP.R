# This script creates Figure 6 in the paper (Section 4)

library(colorspace)
library(ggplot2)


c = c(rgb(1,0.7,0.7), rgb(0.7,1,0.7), rgb(0.7,0.7,1), rgb(1,1,0.7), rgb(1,0.7,1), rgb(0.7,1,1))
mixcolor(0.7, RGB(255,255,0), RGB(255,255,255))

loadings = c(0.3, 0.6, 0.5)
x = matrix(c(1,-1,2,-2,1.5,3,2.3,4,-1,-6.5,8,0,0.5,1,2,4,-2,2), ncol=3,byrow = TRUE)
PC = x %*% loadings

table = data.frame("PC1" = PC, "Mean Prediction" = c(0.5, 0.65, 0.4, 0.3, 0.55, 0.6))
ggplot(data = table, aes(x = PC1, y = Mean.Prediction)) + geom_point(size = 5, colour = c) + 
  theme_bw() + ylab("Mean Prediction") + xlab(expression("PC1 = 0.3x"[1]*" + 0.6x"[2]*" + 0.5x"[3]))
