### tests methods argument of lmrob.control

library(robustbase)

data(stackloss)

## S
set.seed(0)
m0 <- lmrob(stack.loss ~ ., data = stackloss, method = "S",
                    compute.outlier.stats = "S")
m0

write.csv   (m0$coefficients, file = "coeff_Smethod_base.csv", row.names = FALSE )
write.csv(m0$residuals, file = "res_Smethod_base.csv", row.names = FALSE )

## MM
set.seed(0)
m1 <- lmrob(stack.loss ~ ., data = stackloss, method = "MM",
                    compute.outlier.stats = "S")
m1

write.csv(m1$coefficients, file = "coeff_MMmethod_base.csv", row.names = FALSE )
write.csv(m1$residuals, file = "res_MMmethod_base.csv", row.names = FALSE )

## SMD
set.seed(0)
m4 <- lmrob(stack.loss ~ ., data = stackloss, method = "SMD",
                    compute.outlier.stats = "S")
m4

write.csv(m4$coefficients, file = "coeff_SMDmethod_base.csv", row.names = FALSE )
write.csv(m4$residuals, file = "res_SMDmethod_base.csv", row.names = FALSE )

## SMDM
set.seed(0)
m5 <- lmrob(stack.loss ~ ., data = stackloss, method = "SMDM",
                    compute.outlier.stats = "S")
m5

write.csv(m5$coefficients, file = "coeff_SMDMmethod_base.csv", row.names = FALSE )
write.csv(m5$residuals, file = "res_SMDMmethod_base.csv", row.names = FALSE )
