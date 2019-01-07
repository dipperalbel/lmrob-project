library(robustbase)


source(system.file("test-tools-1.R",  package = "Matrix", mustWork=TRUE))
DNase1 <- DNase[ DNase$Run == 1, ]
form <- density ~ Asym/(1 + exp(( xmid -log(conc) )/scal))
pnms <- c("Asym", "xmid", "scal")
psNms <- c(pnms, "sigma")
setNames. <- function(x, nm) setNames(rep_len(x, length(nm)), nm)

write.csv(DNase1, "Function=NLROBInput.csv")
## robust
set.seed(2345) # for reproducibility
Rfit.MM.S.bisquare <-
  nlrob(form, data = DNase1,
           lower = setNames.(0, pnms), upper = c(1),
           method="MM")
write.csv(Rfit.MM.S.bisquare$residuals, "Function=NLROBMM.Data=residuals.csv")
write.csv(Rfit.MM.S.bisquare$coefficients, "Function=NLROBMM.Data=coefficients.csv")

set.seed(2345) # for reproducibility
Rfit.tau.bisquare <- nlrob(form, data = DNase1,
        lower = setNames.(0, pnms), upper = c(1),
        method="tau")

set.seed(2345) # for reproducibility
write.csv(Rfit.tau.bisquare$residuals, "Function=NLROBtau.Data=residuals.csv")
write.csv(Rfit.tau.bisquare$coefficients, "Function=NLROBtau.Data=coefficients.csv")

set.seed(2345) # for reproducibility
Rfit.CM <- nlrob(form, data = DNase1,
                 lower = setNames.(0, pnms), upper = c(1),
                 method="CM")
write.csv(Rfit.CM$residuals, "Function=NLROBCM.Data=residuals.csv")
write.csv(Rfit.CM$coefficients, "Function=NLROBCM.Data=coefficients.csv")

set.seed(2345) # for reproducibility
Rfit.mtl <- nlrob(form, data = DNase1,
                  lower = setNames.(0, pnms), upper = c(1),
                  method="mtl")
write.csv(Rfit.mtl$residuals, "Function=NLROBmtl.Data=residuals.csv")
write.csv(Rfit.mtl$coefficients, "Function=NLROBmtl.Data=coefficients.csv")

Cfit <- nlrob(form, data = DNase1, start= setNames.(1, pnms),
              method="M")
    
write.csv(Cfit$residuals, "Function=NLROBM.Data=residuals.csv")

write.csv(Cfit$coefficients, "Function=NLROBM.Data=coefficients.csv")

