library(robustbase)
nlrob.MM  <- robustbase:::nlrob.MM
nlrob.tau <- robustbase:::nlrob.tau
nlrob.CM  <- robustbase:::nlrob.CM
nlrob.mtl <- robustbase:::nlrob.mtl


source(system.file("test-tools-1.R",  package = "Matrix", mustWork=TRUE))
DNase1 <- DNase[ DNase$Run == 1, ]
form <- density ~ Asym/(1 + exp(( xmid -log(conc) )/scal))
pnms <- c("Asym", "xmid", "scal")
psNms <- c(pnms, "sigma")
setNames. <- function(x, nm) setNames(rep_len(x, length(nm)), nm)

## robust
set.seed(2345) # for reproducibility
Rfit.MM.S.bisquare <-
  nlrob.MM(form, data = DNase1,
           lower = setNames.(0, pnms), upper = c(1),
           ctrl = nlrob.control("MM", optim.control = list(trace = 1),
                                optArgs = list(trace = FALSE)))
write.csv(DNase1, "Submodule=NLROBMM.Data=Input.csv")
write.csv(Rfit.MM.S.bisquare$residuals, "Submodule=NLROBMM.Data=residuals.csv")
write.csv(Rfit.MM.S.bisquare$coefficients, "Submodule=NLROBMM.Data=coefficients.csv")

set.seed(2345) # for reproducibility
Rfit.tau.bisquare <-
         nlrob.tau(form, data = DNase1,
                   lower = setNames.(0, pnms), upper = c(1),
                   ctrl = nlrob.control("tau", optim.control = list(trace = 1),
                                        optArgs = list(trace = FALSE)))

set.seed(2345) # for reproducibility
write.csv(DNase1, "Submodule=NLROBTAU.Data=Input.csv")
write.csv(Rfit.tau.bisquare$residuals, "Submodule=NLROBTAU.Data=residuals.csv")
write.csv(Rfit.tau.bisquare$coefficients, "Submodule=NLROBTAU.Data=coefficients.csv")

set.seed(2345) # for reproducibility
Rfit.CM <- nlrob.CM(form, data = DNase1,
                    lower = setNames.(0, pnms), upper = c(1),
                    ctrl = nlrob.control("CM", optim.control = list(trace = 1),
                                         optArgs = list(trace = FALSE)))
write.csv(DNase1, "Submodule=NLROBCM.Data=Input.csv")
write.csv(Rfit.CM$residuals, "Submodule=NLROBCM.Data=residuals.csv")
write.csv(Rfit.CM$coefficients, "Submodule=NLROBCM.Data=coefficients.csv")

set.seed(2345) # for reproducibility
Rfit.mtl <- nlrob.mtl(form, data = DNase1,
                      lower = setNames.(0, pnms), upper = c(1),
                      ctrl = nlrob.control("mtl", optim.control = list(trace = 1),
                                           optArgs = list(trace = FALSE)))
write.csv(DNase1, "Submodule=NLROBMTL.Data=Input.csv")
write.csv(Rfit.mtl$residuals, "Submodule=NLROBMTL.Data=residuals.csv")
write.csv(Rfit.mtl$coefficients, "Submodule=NLROBMTL.Data=coefficients.csv")

Cfit <- nls(form, data=DNase1, start= c(Asym = 1, scal = 0.1, xmid = 0.4),
            control = nls.control(tol = 8e-8, printEval = TRUE))

write.csv(DNase1, "Submodule=NLS.Data=Input.csv")
write.csv(Cfit$m$resid(), "Submodule=NLS.Data=residuals.csv")
write.csv(Cfit$m$getPars(), "Submodule=NLS.Data=coefficients.csv")

    