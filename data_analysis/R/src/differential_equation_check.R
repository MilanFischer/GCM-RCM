# Clear the workspace
rm(list=ls())

# library(Deriv)
# library(Ryacas)
# library(calculus)

#-------------------------------------------------------------------------------

# ET = geff * VPD

ET_hist <- 500
ET_fut <- 520
d_ET <- ET_fut - ET_hist
VPD_hist <- 0.4
VPD_fut <- 0.8
d_VPD <- VPD_fut - VPD_hist
geff_hist <- ET_hist / VPD_hist
geff_fut <- ET_fut / VPD_fut
d_geff <- geff_fut - geff_hist


# Normalized differences
d_ET_over_ET <- d_ET / ET_hist
d_VPD_over_VPD <- d_VPD / VPD_hist
d_geff_over_geff <- d_geff / geff_hist

# First order Taylor-series expansion
# ET = geff * VPD

#------------------
# First order terms
D(expression(geff * VPD), 'geff')
D(expression(geff * VPD), 'VPD')

#-------------------
# Second order terms
D(D(expression(geff * VPD), 'geff'), 'geff')
D(D(expression(geff * VPD), 'VPD'), 'VPD')

# Cross terms
D(D(expression(geff * VPD), 'geff'), 'VPD')

d_ET_est <-
  1 / factorial(1) * (VPD_hist * d_geff + geff_hist * d_VPD) +
  1 / factorial(2) * (0 * (d_geff^2) + 0 * (d_VPD^2) + 2 * 1 * d_geff * d_VPD)

#-------------------------------------------------------------------------------

# geff = ET / VPD

#------------------
# First order terms
D(expression(ET / VPD), 'ET')
D(expression(ET / VPD), 'VPD')

#-------------------
# Second order terms
D(D(expression(ET / VPD), 'ET'), 'ET')
D(D(expression(ET / VPD), 'VPD'), 'VPD')

# Cross terms
D(D(expression(ET / VPD), 'ET'), 'VPD')

d_geff_est <-
  1 / factorial(1) * (1 / VPD_hist * d_ET + -(ET_hist / VPD_hist^2) * d_VPD) +
  1 / factorial(2) * (0 * (d_ET^2) + ET_hist * (2 * VPD_hist) / (VPD_hist^2)^2 * (d_VPD^2) + 2 * -(1 / VPD_hist^2) * d_ET * d_VPD)

#-------------------------------------------------------------------------------