df <- read.delim("C:/Users/lbolwin/OneDrive - IW/Desktop/Diss/Kurse/Stochastic Processes/PS/R/Exercise 20.txt")

trgroup = df[1:10,2]
plgroup = df[11:21,2]

teststat = mean(trgroup) - mean(plgroup)
t.test(trgroup, plgroup, alternative = "greater")

# Idea for using Bootstrap
# Bootstrap the test statistic and use a percentile confidence interval as decision rule

B = 20000
teststatB = rep(0,B)
set.seed(161093)

for (b in 1:B){
  trB = sample(trgroup, size = 10, replace = T)
  plB = sample(plgroup, size = 11, replace = T)
  
  teststatB[b] = mean(trB) - mean(plB)
}

teststatB = as.data.frame(teststatB)

library(ggplot2)

ggplot(teststatB) +
  aes(x = teststatB) +
  geom_density(adjust = 1L, fill = "#112446", alpha = 0.2) +
  geom_vline(xintercept = quantile(teststatB$teststatB, 0.05), linetype="dotted", 
             color = "red", size=1)+
  geom_vline(xintercept = quantile(teststatB$teststatB, 0.95), linetype="dotted", 
             color = "red", size=1)+
  theme_minimal()

quantile(teststatB$teststatB, 1)
