## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 4.5
)
Sys.setenv(TZ = "UTC")

## ----install-cran, eval=FALSE-------------------------------------------------
# install.packages("INLAcircular")

## ----install-inla, eval=FALSE-------------------------------------------------
# install.packages(
#   "INLA",
#   repos = c(
#     getOption("repos"),
#     INLA = "https://inla.r-inla-download.org/R/stable"
#   )
# )

## ----install-package, eval=FALSE----------------------------------------------
# install.packages("remotes")
# remotes::install_github(
#   "XiangYEstats/INLAcircular"
# )

## ----load---------------------------------------------------------------------
library(INLAcircular)

## ----density-comparison, fig.cap="Three circular density families with mean direction zero."----
theta <- seq(-pi, pi, length.out = 361)
densities <- cbind(
  "von Mises" = dvm(theta, mu = 0, kappa = 2),
  "cardioid" = dcardioid(theta, mu = 0, kappa = 0.4),
  "wrapped Cauchy" =
    dwrappedcauchy(theta, mu = 0, kappa = 0.7)
)

matplot(
  theta,
  densities,
  type = "l",
  lty = 1,
  xlab = "Angle (radians)",
  ylab = "Density"
)
legend(
  "topright",
  legend = colnames(densities),
  col = seq_len(ncol(densities)),
  lty = 1,
  bty = "n"
)

## ----von-mises-api------------------------------------------------------------
set.seed(1)
rvm(5, mu = 0, kappa = 3)
pvm(pi / 2, mu = 0, kappa = 2)
qvm(0.95, mu = 0, kappa = 2)

## ----lavm-api-----------------------------------------------------------------
set.seed(2)
rlavm(5, eta = 1.5, kappa = 3)
dlavm(c(-1, 0, 1), eta = 1.5, kappa = 3)
plavm(pi / 2, eta = 1.5, kappa = 3)
qlavm(0.5, eta = 1.5, kappa = 3)

## ----other-densities----------------------------------------------------------
x <- c(-pi / 2, 0, pi / 2)
dcardioid(x, mu = 0, kappa = c(0, 0.25, 0.5))
dwrappedcauchy(x, mu = 0, kappa = 0.7, log = TRUE)

## ----pc-vm--------------------------------------------------------------------
probability <- c(0.1, 0.5, 0.9)
q0 <- qpc.vm0(probability, lambda = 2)
cbind(
  probability = probability,
  quantile = q0,
  recovered = ppc.vm0(q0, lambda = 2)
)

dpc.vm0(c(0, 0.5, 2), lambda = 2)

## ----pc-vminf-----------------------------------------------------------------
lambda <- 2
atom_vm <- exp(-lambda)
c(computed = ppc.vminf(0, lambda), expected = atom_vm)
dpc.vminf(c(0, 0.5, 2), lambda)
qpc.vminf(c(atom_vm / 2, (1 + atom_vm) / 2), lambda)

## ----pc-card0-----------------------------------------------------------------
probability <- c(0.1, 0.5, 0.9)
quantile <- qpc.card0(probability, lambda = 2)
cbind(
  probability = probability,
  quantile = quantile,
  recovered = ppc.card0(quantile, lambda = 2)
)
dpc.card0(c(0, 0.25, 0.5), lambda = 2)

## ----pc-card------------------------------------------------------------------
lambda <- 2
atom_card <- exp(-lambda * sqrt(log(2)))
c(computed = ppc.card(0, lambda), expected = atom_card)
qpc.card(c(atom_card / 2, (1 + atom_card) / 2), lambda)
dpc.card(c(0, 0.25, 0.49), lambda)

## ----pc-wc--------------------------------------------------------------------
q_wc <- qpc.wc(c(0.1, 0.5, 0.9), lambda = 2)
cbind(quantile = q_wc, recovered = ppc.wc(q_wc, lambda = 2))

set.seed(3)
rpc.wc(5, lambda = 2)
dpc.wc(c(0, 0.5, 0.9), lambda = 2)

## ----pc-log-vector------------------------------------------------------------
dpc.card0(c(0.1, 0.2, 0.3), lambda = c(1, 2, 3), log = TRUE)
qpc.wc(log(c(0.1, 0.5, 0.9)), lambda = 2, log.p = TRUE)

## ----direct-inla, eval=FALSE--------------------------------------------------
# lavm_setting <- list(
#   link = "inverse.tangent",
#   hyper = list(
#     kappa = list(
#       initial = log(50),
#       prior = "pc.vminf",
#       param = c(0.5, 0.5),
#       fixed = FALSE
#     )
#   )
# )
# 
# direct_data <- list(
#   angle_lavm = INLA::inla.mdata(angle),
#   z = z
# )
# 
# fit_direct <- INLA::inla(
#   angle_lavm ~ z,
#   family = "cloglike",
#   data = direct_data,
#   control.family = list(
#     cloglike = lavm.cloglike(lavm_setting)
#   ),
#   control.inla = list(
#     cmin = 0,
#     compute.initial.values = TRUE
#   ),
#   control.compute = list(cpo = TRUE, dic = TRUE, waic = TRUE),
#   verbose = FALSE
# )

## ----direct-inla-transform, eval=FALSE----------------------------------------
# theta_name <- grep(
#   "^Theta1 for",
#   names(fit_direct$marginals.hyperpar),
#   value = TRUE
# )[1]
# kappa_marginal <- INLA::inla.tmarginal(
#   exp,
#   fit_direct$marginals.hyperpar[[theta_name]]
# )
# INLA::inla.zmarginal(kappa_marginal)

## ----convenient-inla, eval=FALSE----------------------------------------------
# library(INLA)
# library(INLAcircular)
# 
# fit_short <- inla(
#   angle ~ z,
#   family = "lavm",
#   data = data.frame(angle = angle, z = z),
#   control.family = lavm_setting
# )

## ----model-specification, eval=FALSE------------------------------------------
# angle_model <- likelihood(
#   angle ~
#     intercept(name = "alpha0", mean = 0, sd = 2) +
#     covariate(z, name = "alpha_z", mean = 0, sd = 1) +
#     f(
#       hour_effect,
#       model = "rw1",
#       cyclic = TRUE,
#       constr = TRUE,
#       scale.model = TRUE,
#       hyper = list(
#         prec = list(
#           prior = "pc.prec",
#           param = c(1, 0.01)
#         )
#       )
#     ),
#   family = "lavm",
#   family.setting = lavm_setting
# )

## ----latent-index-------------------------------------------------------------
hour_index <- index(
  var = "hour_effect",
  data.id = c(1, 2, 2, 4),
  process.id = 1:4
)

hour_index

## ----multiple-indices, eval=FALSE---------------------------------------------
# latent_indices <- list(
#   index("hour_effect", data.id = data$hour, process.id = 0:23),
#   index("region_effect", data.id = data$region, process.id = region_ids)
# )

## ----inlacc-fit, eval=FALSE---------------------------------------------------
# fit <- inlacc(
#   model = angle_model,
#   data = data.frame(
#     angle = angle,
#     z = z,
#     hour_effect = hour
#   ),
#   latent.index = index(
#     "hour_effect",
#     data.id = hour,
#     process.id = 0:23
#   ),
#   control.predictor = list(compute = TRUE),
#   metrics = TRUE
# )

## ----joint-inlacc, eval=FALSE-------------------------------------------------
# joint_model <- list(
#   likelihood(
#     angle ~
#       intercept(name = "alpha0", mean = 0, sd = 2) +
#       covariate(z, name = "alpha_z", mean = 0, sd = 1),
#     family = "lavm",
#     family.setting = lavm_setting
#   ),
#   likelihood(
#     y ~
#       intercept(name = "beta0", mean = 0, sd = 5) +
#       covariate(z, name = "beta_z", mean = 0, sd = 2) +
#       covariate(
#         angle,
#         name = "beta_angle_predictor",
#         mean = 0,
#         sd = 1,
#         initial = 0,
#         fixed = FALSE
#       ),
#     family = "gaussian",
#     family.setting = list(
#       hyper = list(
#         prec = list(prior = "pc.prec", param = c(1, 0.01))
#       )
#     )
#   )
# )
# 
# joint_fit <- inlacc(
#   model = joint_model,
#   data = data.frame(angle = angle, y = y, z = z),
#   control.predictor = list(compute = TRUE),
#   metrics = TRUE
# )

## ----graphpcor-install, eval=FALSE--------------------------------------------
# install.packages(c("graphpcor", "INLAtools"))

## ----graphpcor-model, eval=FALSE----------------------------------------------
# multivariate_model <- list(
#   likelihood(
#     angle ~ intercept(name = "alpha0") +
#       f(
#         i,
#         model = "iidkd_LKJ",
#         pc.prior.u = 1,
#         pc.prior.alpha = 0.5
#       ),
#     family = "lavm",
#     family.setting = lavm_setting
#   ),
#   likelihood(
#     y ~ intercept(name = "beta0") +
#       f(
#         i,
#         model = "iidkd_LKJ",
#         pc.prior.u = 2,
#         pc.prior.alpha = 0.5
#       ),
#     family = "gaussian"
#   )
# )
# 
# fit_lkj <- inlacc(
#   multivariate_model,
#   data = data.frame(angle = angle, y = y),
#   LKJ.eta = 5,
#   control.predictor = list(compute = TRUE)
# )

## ----graphpcor-index, eval=FALSE----------------------------------------------
# lkj_index <- index(
#   var = "i",
#   data.id = replicate_id,
#   process.id = "likelihood"
# )
# 
# fit_lkj <- inlacc(
#   multivariate_model,
#   data = data,
#   latent.index = lkj_index,
#   LKJ.eta = 5
# )

## ----graphpcor-reconstruct, eval=FALSE----------------------------------------
# p <- 2L
# theta_i <- fit_lkj$lkj_modes[["i"]]
# 
# sd_mode <- exp(theta_i[seq_len(p)])
# correlation_mode <- graphpcor::basecor(
#   theta_i[-seq_len(p)],
#   p = p
# )$base
# 
# sd_mode
# correlation_mode

## ----graphpcor-reconstruct-fixed, eval=FALSE----------------------------------
# process_i <- fit_lkj$inlacc_meta$LKJ[["i"]]
# u <- process_i$pc.prior.u
# alpha <- process_i$pc.prior.alpha
# fixed_scale <- is.na(alpha) | alpha <= 0 | alpha >= 1
# n_free <- sum(!fixed_scale)
# 
# theta_i <- fit_lkj$lkj_modes[["i"]]
# log_sd_mode <- log(u)
# log_sd_mode[!fixed_scale] <- theta_i[seq_len(n_free)]
# 
# n_cor <- p * (p - 1L) / 2L
# cor_index <- n_free + seq_len(n_cor)
# correlation_mode <- graphpcor::basecor(
#   theta_i[cor_index],
#   p = p
# )$base
# 
# sd_mode <- exp(log_sd_mode)

## ----summaries, eval=FALSE----------------------------------------------------
# fit_summary <- summary(fit, decimal = 3)
# print(fit_summary, metrics = TRUE)

## ----bessel-------------------------------------------------------------------
bessel_i(c(0, 1, 5), nu = 0)
bessel_i(1, nu = 0.5)
bessel_i(1000, nu = 0, expon.scaled = TRUE)

## ----wind-data----------------------------------------------------------------
data("wind_newyork", package = "INLAcircular")
dim(wind_newyork)
utils::head(names(wind_newyork), 4)

## ----session-info-------------------------------------------------------------
R.version.string
packageVersion("INLAcircular")

