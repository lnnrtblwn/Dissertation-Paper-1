# Testen. Was bringen die Covariates in den Anwendungen überhaupt?
# irgendwann dann machen

synth(X1 = matrix(c(y_before[1:(T0 %/% 2),1], dummies[1,])),  
      X0 = rbind(y_before[1:(T0 %/% 2),-1], t(dummies[-1,])), 
      Z1 = matrix(y_before[,1]),  
      Z0 = y_before[,-1])

# was machen ADH?

data("basque")

basque[85:89, 1:4]

dataprep.out <- dataprep(
  foo=basque,
  predictors=c("school.illit","school.prim","school.med","school.high", "school.post.high", "invest"),
  predictors.op="mean",
  time.predictors.prior=1964:1969,
  special.predictors=list(
    list("gdpcap", 1960:1969 , "mean"),
    list("sec.agriculture", seq(1961, 1969, 2), "mean"),
    list("sec.energy", seq(1961, 1969, 2), "mean"),
    list("sec.industry", seq(1961, 1969, 2), "mean"),
    list("sec.construction", seq(1961, 1969, 2), "mean"),
    list("sec.services.venta", seq(1961, 1969, 2), "mean"),
    list("sec.services.nonventa", seq(1961, 1969, 2), "mean"),
    list("popdens", 1969, "mean")),
  dependent="gdpcap",
  unit.variable="regionno",
  unit.names.variable="regionname",
  time.variable="year",
  treatment.identifier=17,
  controls.identifier=c(2:16,18),
  time.optimize.ssr=1960:1969,
  time.plot=1955:1997)

dataprep.out$X0
dataprep.out$X1
dataprep.out$Z0
dataprep.out$Z1

quiet <- function(x) { 
  sink(tempfile()) 
  on.exit(sink()) 
  invisible(force(x)) 
} 
