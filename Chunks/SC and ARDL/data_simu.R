
library(reshape2)

### Dataset Simulation

set.seed(1)

year <- 1970:2022
cutoff = 2000
nElem = 5


# with covariance
library(MASS)
#Sigma <- matrix(c(10,3,3,3,3,5,3,3,3,3,8,3,3,3,3,15)/10,ncol=nvars)
Sigma <- rWishart(n=1, df=nElem, Sigma=diag(nvElem))[,,1]
Sigma
smdata=mvrnorm(n = length(year), rep(0, nElem), Sigma)
smcdata = cumsum(data.frame(smdata))
scdf=data.frame(year,smcdata)
scdfy = scdf %>% 
  mutate(Y= 0.2*X1 + 0.5 * X2 + 0.1 * X3 + 0.15 * X4 + 0.05*X5)
summary(scdfy)

Sigma
cov(scdf)





#Variable Gen

df=data.frame(rep(year,nElem))

set.seed(1)
Sigma <- rWishart(n=1, df=nElem, Sigma=diag(nElem))[,,1]
V1=data.frame(mvrnorm(n = length(year), rep(0, nElem), Sigma))
colnames(V1) = c("A","B","C","D","E")
df=data.frame(df,reshape2::melt(V1))

set.seed(2)
Sigma <- rWishart(n=1, df=nElem, Sigma=diag(nElem))[,,1]
V2=data.frame(mvrnorm(n = length(year), rep(0, nElem), Sigma))
colnames(V2) = c("A","B","C","D","E")
df=data.frame(df,reshape2::melt(V2)[2])

set.seed(3)
Sigma <- rWishart(n=1, df=nElem, Sigma=diag(nElem))[,,1]
V3=data.frame(mvrnorm(n = length(year), rep(0, nElem), Sigma))
colnames(V3) = c("A","B","C","D","E")
df=data.frame(df,reshape2::melt(V3)[2])

set.seed(4)
Sigma <- rWishart(n=1, df=nElem, Sigma=diag(nElem))[,,1]
V4=data.frame(mvrnorm(n = length(year), rep(0, nElem), Sigma))
colnames(V4) = c("A","B","C","D","E")
df=data.frame(df,reshape2::melt(V4)[2])

Y = 0.4*V1 + 0.35*V2 + +0.2*V3+0.05*V4 

colnames(Y) = c("A","B","C","D","E")

df=data.frame(df,reshape2::melt(Y)[2])

colnames(df) = c("Year","Elem","X1","X2","X3","X4","Y")



df = df %>% 
  mutate(treat = )

### Treatment Faktor
treat_fak = rep(c(0,1,0),times=c(31,22,212))
scdft = df %>% 
  mutate(Y_t =  Y+(5*treat_fak))

ggplot(reshape2::melt(scdft,id="year")) +
  aes(x = year, y = value, colour = variable) +
  geom_line(size = 0.75) +
  scale_color_hue(direction = 1) +
  geom_vline(xintercept = 2000, linetype="dashed", size=.5)+
  theme_minimal()





