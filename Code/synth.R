library("Synth")
library(dplyr)
data("basque")

basque[85:89, 1:4]

dataprep.out <- dataprep(
  foo = basque,
  predictors = c("school.illit", "school.prim", "school.med",
                 "school.high", "school.post.high", "invest"),
  predictors.op = "mean",
  time.predictors.prior = 1964:1969,
  special.predictors = list(
    list("gdpcap", 1956:1969 , "mean"),
    list("sec.agriculture", seq(1961, 1969, 2), "mean"),
    list("sec.energy", seq(1961, 1969, 2), "mean"),
    list("sec.industry", seq(1961, 1969, 2), "mean"),
    list("sec.construction", seq(1961, 1969, 2), "mean"),
    list("sec.services.venta", seq(1961, 1969, 2), "mean"),
    list("sec.services.nonventa", seq(1961, 1969, 2), "mean"),
    list("popdens", 1969, "mean")),
  dependent = "gdpcap",
  unit.variable = "regionno",
  unit.names.variable = "regionname",
  time.variable = "year",
  treatment.identifier = 17,
  # treatment.identifier = 10, # placebo
  controls.identifier = c(2:16, 18),
  # controls.identifier = c(2:9, 11:16, 18), # placebo
  time.optimize.ssr = 1960:1969,
  time.plot = 1955:1997)

dataprep.out$X0
dataprep.out$X1

test = basque %>% 
  filter(regionname == "Basque Country (Pais Vasco)",
         year %in% c(1964:1969))
mean(test$invest)
rm(test)

dataprep.out$Z1

# some data manipulation

dataprep.out$X1["school.high",] <- dataprep.out$X1["school.high",] + 
  dataprep.out$X1["school.post.high",]

dataprep.out$X1 <- as.matrix(dataprep.out$X1[
  + -which(rownames(dataprep.out$X1) == "school.post.high"),])

dataprep.out$X0["school.high",] <- dataprep.out$X0["school.high",] +
  + dataprep.out$X0["school.post.high",]

dataprep.out$X0 <- dataprep.out$X0[
  + -which(rownames(dataprep.out$X0) == "school.post.high"),]

lowest <- which(rownames(dataprep.out$X0) == "school.illit")
highest <- which(rownames(dataprep.out$X0) == "school.high")
dataprep.out$X1[lowest:highest,] <-
  + (100 * dataprep.out$X1[lowest:highest,]) /
  + sum(dataprep.out$X1[lowest:highest,])
dataprep.out$X0[lowest:highest,] <-
  100 * scale(dataprep.out$X0[lowest:highest,], center = FALSE,
              scale = colSums(dataprep.out$X0[lowest:highest,]))

# Running synth

synth.out <- synth(data.prep.obj = dataprep.out, method = "BFGS")

gaps <- dataprep.out$Y1plot - (dataprep.out$Y0plot %*% synth.out$solution.w)
gaps[1:3, 1]

synth.tables <- synth.tab(dataprep.res = dataprep.out,
                          synth.res = synth.out)

synth.tables$tab.pred[1:5, ]
synth.tables$tab.v
synth.tables$tab.w

path.plot(synth.res = synth.out, dataprep.res = dataprep.out,
          Ylab = "real per-capita GDP (1986 USD, thousand)", Xlab = "year",
          Ylim = c(0, 12), Legend = c("Basque country", "synthetic Basque country"), Legend.position = "bottomright")

gaps.plot(synth.res = synth.out, dataprep.res = dataprep.out,
          Ylab = "gap in real per-capita GDP (1986 USD, thousand)", Xlab = "year",
          Ylim = c(-1.5, 1.5), Main = NA)

