# Unit conversions
W_to_mm <- 0.408*3600*24/10^6

# Psychrometric constant
gamma <- 0.063

# Clausius-Clapeyron parameters (Katul et al. 2012)
a <- 0.6111
b <- 17.5
c <- 240.97

# # Clausius-Clapeyron parameters (Savage 2009)
# a <- 0.6108
# b <- 17.2694
# c <- 237.3

# b*c* e_sat/(c + Tair)^2
# b*c* (a*exp(b*Tair/(Tair+c)))/(c + Tair)^2

# Priestley-Taylor coefficient
alpha_PT_1.26 <- 1.26

# The list of models
models <- c("CMIP5", "CMIP6", "RCMs")

# The list of periods
periods <- c('1981_2005', '2076_2100')

days_in_month <- c(31,28.25,31,30,31,30,31,31,30,31,30,31)

##############################################
# The logic of this needs to be double-checked

# Function for annual sum
annSums <- function(x){
  
  if(nrow(x) != 12){
    stop('There must be excatly 12 lines each representing one month\n')
  }
  x[1,] <- x[1,]*31
  x[2,] <-x[2,]*28.25
  x[3,] <-x[3,]*31
  x[4,] <-x[4,]*30
  x[5,] <-x[5,]*31
  x[6,] <-x[6,]*30
  x[7,] <-x[7,]*31
  x[8,] <-x[8,]*31
  x[9,] <-x[9,]*30
  x[10,] <-x[10,]*31
  x[11,] <-x[11,]*30
  x[12,] <-x[12,]*31
  y <- 12* colMeans(x, na.rm=TRUE)/mean(days_in_month)
  
  return(y)
}

# Function for annual sum
annMeans <- function(x){
  
  if(nrow(x) != 12){
    stop('There must be excatly 12 lines each representing one month\n')
  }
  x[1,] <- x[1,]*31
  x[2,] <-x[2,]*28.25
  x[3,] <-x[3,]*31
  x[4,] <-x[4,]*30
  x[5,] <-x[5,]*31
  x[6,] <-x[6,]*30
  x[7,] <-x[7,]*31
  x[8,] <-x[8,]*31
  x[9,] <-x[9,]*30
  x[10,] <-x[10,]*31
  x[11,] <-x[11,]*30
  x[12,] <-x[12,]*31
  y <- colMeans(x, na.rm=TRUE)/mean(days_in_month)
  
  return(y)
}

######################
# Temporal sub-setting
CMIP5_1981_2005_data <- CMIP5_data[which(CMIP5_data$PERIOD=='1981_2005'),]
CMIP5_2076_2100_data <- CMIP5_data[which(CMIP5_data$PERIOD=='2076_2100'),]

CMIP6_1981_2005_data <- CMIP6_data[which(CMIP6_data$PERIOD=='1981_2005'),]
CMIP6_2076_2100_data <- CMIP6_data[which(CMIP6_data$PERIOD=='2076_2100'),]

RCMs_1981_2005_data <- RCMs_data[which(RCMs_data$PERIOD=='1981_2005'),]
RCMs_2076_2100_data <- RCMs_data[which(RCMs_data$PERIOD=='2076_2100'),]

# Select the column with variables
CMIP5_1981_2005 <- 2:(ncol(CMIP5_1981_2005_data)-2)
CMIP5_2076_2100 <- 2:(ncol(CMIP5_2076_2100_data)-2)
CMIP6_1981_2005 <- 2:(ncol(CMIP6_1981_2005_data)-2)
CMIP6_2076_2100 <- 2:(ncol(CMIP6_2076_2100_data)-2)
RCMs_1981_2005 <- 2:(ncol(RCMs_1981_2005_data)-2)
RCMs_2076_2100 <- 2:(ncol(RCMs_2076_2100_data)-2)

############################
# Fluxes and state variables

# Define the function
assign_var <- function(model, var_name, period){
  x_0 <- get(paste0(model, "_", period, "_data"))
  
  if(var_name=="Rn"){
    SW_net <- x_0[which(x_0$VAR=="SW_net"),]
    LW_net <- x_0[which(x_0$VAR=="LW_net"),]
    data_range <- get(paste0(model, "_", period))
    
    # Correct the sign
    LW_net[,data_range] <- (-1)*LW_net[,data_range]
    
    x <- SW_net
    x$VAR <- 'Rn'
    x[,data_range] <- SW_net[,data_range] + LW_net[,data_range]
  }else if(var_name=="LW_net"){
    
    x <- x_0[which(x_0$VAR==var_name),]
    
    # Correct the sign
    data_range <- get(paste0(model, "_", period))
    x[,data_range] <- (-1)*x[,data_range]
    
  }else if(var_name=="G"){
    SW_net <- x_0[which(x_0$VAR=="SW_net"),]
    LW_net <- x_0[which(x_0$VAR=="LW_net"),]
    H <- x_0[which(x_0$VAR=="H"),]
    LE <- x_0[which(x_0$VAR=="LE"),]
    data_range <- get(paste0(model, "_", period))
    
    # Correct the sign
    LW_net[,data_range] <- (-1)*LW_net[,data_range]
    
    x <- SW_net
    x$VAR <- 'G'
    x[,data_range] <- SW_net[,data_range] + LW_net[,data_range] - H[,data_range] - LE[,data_range]
  }else if(var_name=="e_sat"){
    Ta <- x_0[which(x_0$VAR=="Ta"),]
    x <- Ta
    x$VAR <- 'e_sat'
    data_range <- get(paste0(model, "_", period))
    x[,data_range] <- a*exp(b*Ta[,data_range]/(Ta[,data_range]+c))
  }else if(var_name=="SVP"){
    Ta <- x_0[which(x_0$VAR=="Ta"),]
    e_sat <- Ta
    e_sat$VAR <- 'e_sat'
    x <- Ta
    x$VAR <- 'SVP'
    data_range <- get(paste0(model, "_", period))
    e_sat[,data_range] <- a*exp(b*Ta[,data_range]/(Ta[,data_range]+c))
    x[,data_range] <-  b*c*e_sat[,data_range]/(c + Ta[,data_range])^2
  }else{
    x <- x_0[which(x_0$VAR==var_name),]
  }
  
  assign(x = paste0(var_name, "_", model, "_", period) , value = x, envir = parent.frame())
}


##############################
# Evapotranspiration function

# Define the function
assign_ET_var <- function(model, var_name, period){
  x_0 <- get(paste0(model, "_", period, "_data"))
  
  n_days <- matrix(rep(days_in_month, ncol(x_0)-4), nrow=12, ncol=ncol(x_0)-4)
  
  if(var_name=="ET"){
    LE <- x_0[which(x_0$VAR=="LE"),]
    x <- LE
    x$VAR <- 'ET'
    data_range <- get(paste0(model, "_", period))
    
    x[,data_range] <-  x[,data_range]*n_days*W_to_mm
  }else if (var_name=="ET_eq"){
    SVP <- get(paste0("SVP_", model, "_", period))
    H <- get(paste0("H_", model, "_", period))
    LE <- get(paste0("LE_", model, "_", period))
    
    data_range <- get(paste0(model, "_", period))
    x <- LE
    x$VAR <- 'Equilibrium_ET'
    
    x[,data_range] <-  SVP[,data_range]*(H[,data_range]+LE[,data_range])/(SVP[,data_range] + gamma)*n_days*W_to_mm
    
    # Replace negative values by zero
    x[,data_range][x[,data_range] < 0] <- 0
  }else if (var_name=="PET_PT_1.26"){
    ET_eq <- get(paste0("ET_eq_", model, "_", period))
    
    data_range <- get(paste0(model, "_", period))
    x <- ET_eq
    x$VAR <- 'PET_PT_alpha_1.26'
    
    x[,data_range] <- ET_eq[,data_range]*alpha_PT_1.26
  }else if (var_name=="PET_PT"){
    ET_eq_baseline <- get(paste0("ET_eq_", model, "_1981_2005"))
    ET_baseline <- get(paste0("ET_", model, "_1981_2005"))
    
    data_range_baseline <- get(paste0(model, "_1981_2005"))
    
    alpha_PT <- ET_eq_baseline
    alpha_PT$VAR <- 'alpha_PT'
    alpha_PT[,data_range_baseline] <- ET_baseline[,data_range_baseline]/ET_eq_baseline[,data_range_baseline]
    
    alpha_PT <- alpha_PT[alpha_PT$MONTH %in% c('IV','V','VI','VII','VIII'),]
    max_alpha_PT <- apply(alpha_PT[,data_range_baseline], 2, max)
    
    alpha_PT <- alpha_PT[1,]
    alpha_PT[,data_range_baseline] <- max_alpha_PT
    alpha_PT$MONTH <- 'IV-VIII'
    
    assign(x = paste0("alpha_PT_", model, "_1981_2005") , value = alpha_PT, envir = parent.frame())
    
    # Now for any period
    ET_eq <- get(paste0("ET_eq_", model, "_", period))
    x <- get(paste0("ET_eq_", model, "_", period))
    x$VAR <- 'PET_PT_alpha_dynamic'
    
    data_range <- get(paste0(model, "_", period))
    
    for(i in data_range){
      x[,i] <- ET_eq[,i]*max(alpha_PT[i], alpha_PT_1.26)
    }
  }else if(var_name=="PET_PT_Zhou"){
    # Following Zhou et al. (2023) "https://doi.org/10.1038/s41558-023-01659-8" This indicated that the estimated PET from the PET algorithms in equations (1)–(4) is probably underestimated. To overcome this problem, we set monthly PET to be the evapotranspiration if the calculated PET is less than evapotranspiration.
    ET_eq <- get(paste0("ET_eq_", model, "_", period))
    ET <- get(paste0("ET_", model, "_", period))
    
    data_range<- get(paste0(model, "_", period))
    
    x <- ET_eq
    x$VAR <- var_name
    x[,data_range] <- alpha_PT_1.26*ET_eq[,data_range]
    for(i in data_range){
      x[,i] <- ifelse(x[,i]<ET[,i],  ET[,i], x[,i])
    }
    
    alpha_PT_Zhou <- ET_eq
    alpha_PT_Zhou$VAR <- 'alpha_PT_Zhou'
    alpha_PT_Zhou[,data_range] <- ET[,data_range]/x[,data_range]
    
    assign(x = paste0("alpha_PT_Zhou", model, "_", period) , value = alpha_PT_Zhou, envir = parent.frame())
  }else if(var_name=="PET_PT_Morton"){
    H <- get(paste0("H_", model, "_", period))
    LE <- get(paste0("LE_", model, "_", period))
    LW_net <- get(paste0("LW_net_", model, "_", period))
    SVP <- get(paste0("SVP_", model, "_", period))
    Rg <- get(paste0("Rg_", model, "_", period))
    
    data_range <- get(paste0(model, "_", period))
    x <- LE
    x$VAR <- 'PET_PT_Morton'
    
    Mm <- LE
   
    Mm[,data_range] <- (-1.37*LW_net[,data_range] - 0.394*Rg[,data_range])
    
    for(i in data_range){
      Mm[,i] <- ifelse(Mm[,i]<0,  0, Mm[,i])
    }

    x[,data_range] <-  alpha_PT_1.26*SVP[,data_range]*(H[,data_range]+LE[,data_range]+Mm[,data_range])/(SVP[,data_range] + gamma)*n_days*W_to_mm
  }else if(var_name=="PET_Makkink_Hansen"){
    SVP <- get(paste0("SVP_", model, "_", period))
    Rg <- get(paste0("Rg_", model, "_", period))
    
    data_range <- get(paste0(model, "_", period))
    x <- Rg
    x$VAR <- 'PET_Makkink_Hansen'
    
    x[,data_range] <-  (0.7*SVP[,data_range]/(SVP[,data_range] + gamma)*Rg[,data_range])*n_days*W_to_mm
  }else if (var_name=="PET_deBruin_2016"){
    SVP <- get(paste0("SVP_", model, "_", period))
    H <- get(paste0("H_", model, "_", period))
    LE <- get(paste0("LE_", model, "_", period))
    
    data_range <- get(paste0(model, "_", period))
    x <- LE
    x$VAR <- 'PET_deBruin_2016'
    
    # Following de Bruin and Holtslag (1982), the correction consists of adding to LE_eq a constant b of approximately 20 W/m2
    Alpha <- 1.0
    Beta <- 20
    x[,data_range] <-  (Alpha*SVP[,data_range]*(H[,data_range]+LE[,data_range])/(SVP[,data_range] + gamma)+Beta)*n_days*W_to_mm
    
    PT_original <- x
    PT_original[,data_range] <- (1.26*SVP[,data_range]*(H[,data_range]+LE[,data_range])/(SVP[,data_range] + gamma))*n_days*W_to_mm

    # Use the greater value
    x[,data_range] <- pmax(PT_original[,data_range], x[,data_range])
    
    # Replace negative values by zero
    x[,data_range][x[,data_range] < 0] <- 0
  }else if (var_name=="PET_deBruin_1979"){
    SVP <- get(paste0("SVP_", model, "_", period))
    H <- get(paste0("H_", model, "_", period))
    LE <- get(paste0("LE_", model, "_", period))
    
    data_range <- get(paste0(model, "_", period))
    x <- LE
    x$VAR <- 'PET_deBruin_1979'
    
    # Following de Bruin and Keijman (1979)
    Alpha <- 1.17
    Beta <- 10
    x[,data_range] <- (Alpha*SVP[,data_range]*(H[,data_range]+LE[,data_range])/(SVP[,data_range] + gamma)+Beta)*n_days*W_to_mm
    
    # Replace negative values by zero
    x[,data_range][x[,data_range] < 0] <- 0
  }
  
  assign(x = paste0(var_name, "_", model, "_", period) , value = x, envir = parent.frame())
}

# Define the function
# Following Zhou et al. (2023) "https://doi.org/10.1038/s41558-023-01659-8" This indicated that the estimated PET from the PET algorithms in equations (1)–(4) is probably underestimated. To overcome this problem, we set monthly PET to be the evapotranspiration if the calculated PET is less than evapotranspiration.
adjust_PET_annualy <- function(model, var_name, period, adjustement_period="baseline"){
  
  if(var_name=="PET_PT"){
    
    if(adjustement_period=="baseline"){
      # Adjustment only for baseline
      ET_baseline <- get(paste0("ET_", model, "_1981_2005"))
      PET_PT_baseline <- get(paste0("PET_PT_", model, "_1981_2005"))
      ET_eq_baseline <- get(paste0("ET_eq_", model, "_1981_2005"))
      alpha_PT <- get(paste0("alpha_PT_", model, "_1981_2005"))
      
      data_range_baseline <- get(paste0(model, "_1981_2005"))
      
      ET_ann <- annSums(ET_baseline[,data_range_baseline])
      PET_PT_ann <- annSums(PET_PT_baseline[,data_range_baseline])
      ET_eq_ann <- annSums(ET_eq_baseline[,data_range_baseline])
      
      tmp <- alpha_PT 
      tmp[,data_range_baseline] <- ET_ann/ET_eq_ann
      
      for(i in data_range_baseline){
        alpha_PT[,i] <- ifelse(tmp[,i]>alpha_PT[,i], tmp[,i], alpha_PT[,i]) # Ensuring that PET is not computed for alpha smaller than 1.26
      }
      
      # Now for any period
      ET_eq <- get(paste0("ET_eq_", model, "_", period))
      
      data_range <- get(paste0(model, "_", period))
      PET_PT <- ET_eq
      PET_PT$VAR <- 'PET_PT_alpha_dynamic'
      
      for(i in data_range){
        PET_PT[,i] <- ET_eq[,i]*max(alpha_PT[i], alpha_PT_1.26)
      }
      
      assign(x = paste0("alpha_PT_", model, "_1981_2005") , value = alpha_PT, envir = parent.frame())
      assign(x = paste0("PET_PT_", model, "_", period) , value = PET_PT, envir = parent.frame())
    }else if(adjustement_period=="both"){
      # Adjustment for baseline and future period
      ET <- get(paste0("ET_", model, "_", period))
      PET_PT <- get(paste0("PET_PT_", model, "_", period))
      ET_eq <- get(paste0("ET_eq_", model, "_", period))
      alpha_PT <- get(paste0("alpha_PT_", model, "_1981_2005"))
      
      data_range <- get(paste0(model, "_", period))
      
      ET_ann <- annSums(ET[,data_range])
      PET_PT_ann <- annSums(PET_PT[,data_range])
      ET_eq_ann <- annSums(ET_eq[,data_range])
      
      tmp <- alpha_PT 
      tmp[,data_range] <- ET_ann/ET_eq_ann
      
      for(i in data_range){
        alpha_PT[,i] <- ifelse(tmp[,i]>alpha_PT[,i], tmp[,i], alpha_PT[,i])
        PET_PT[,i] <- ET_eq[,i]*max(alpha_PT[i], alpha_PT_1.26) # Ensuring that PET is not computed from alpha smaller than 1.26
      }
      
      assign(x = paste0("alpha_PT_", model, "_", period) , value = alpha_PT, envir = parent.frame())
      assign(x = paste0("PET_PT_", model, "_", period) , value = PET_PT, envir = parent.frame())
    }else{
      stop('Adjustement period must be set either to "baseline" or "both"\n')
    }
    
  }
}

# Loop over the models
for(M in models){
  # LE
  assign_var(model = M, var_name = "LE", period = "1981_2005")
  assign_var(model = M, var_name = "LE", period = "2076_2100")
  
  # H
  assign_var(model = M, var_name = "H", period = "1981_2005")
  assign_var(model = M, var_name = "H", period = "2076_2100")
  
  # Rg
  assign_var(model = M, var_name = "Rg", period = "1981_2005")
  assign_var(model = M, var_name = "Rg", period = "2076_2100")
  
  # SW_net
  assign_var(model = M, var_name = "SW_net", period = "1981_2005")
  assign_var(model = M, var_name = "SW_net", period = "2076_2100")
  
  # LW_net
  assign_var(model = M, var_name = "LW_net", period = "1981_2005")
  assign_var(model = M, var_name = "LW_net", period = "2076_2100")
  
  # Rn
  assign_var(model = M, var_name = "Rn", period = "1981_2005")
  assign_var(model = M, var_name = "Rn", period = "2076_2100")
  
  # G
  assign_var(model = M, var_name = "G", period = "1981_2005")
  assign_var(model = M, var_name = "G", period = "2076_2100")
  
  # P
  assign_var(model = M, var_name = "P", period = "1981_2005")
  assign_var(model = M, var_name = "P", period = "2076_2100")
  
  # Ta
  assign_var(model = M, var_name = "Ta", period = "1981_2005")
  assign_var(model = M, var_name = "Ta", period = "2076_2100")
  
  # Saturated vapor pressure
  assign_var(model = M, var_name = "e_sat", period = "1981_2005")
  assign_var(model = M, var_name = "e_sat", period = "2076_2100")
  
  # Slope of the saturated vapor pressure curve
  assign_var(model = M, var_name = "SVP", period = "1981_2005")
  assign_var(model = M, var_name = "SVP", period = "2076_2100")
  
  # Actual evapotranspiration
  assign_ET_var(model = M, var_name = "ET", period = "1981_2005")
  assign_ET_var(model = M, var_name = "ET", period = "2076_2100")
  
  # Equilibrium evaporation
  assign_ET_var(model = M, var_name = "ET_eq", period = "1981_2005")
  assign_ET_var(model = M, var_name = "ET_eq", period = "2076_2100")
  
  # Priestley-Taylor potential evapotranspiration with alpha set to 1.26
  assign_ET_var(model = M, var_name = "PET_PT_1.26", period = "1981_2005")
  assign_ET_var(model = M, var_name = "PET_PT_1.26", period = "2076_2100")
  
  # Priestley-Taylor potential evapotranspiration with dynamic alpha
  assign_ET_var(model = M, var_name = "PET_PT", period = "1981_2005")
  assign_ET_var(model = M, var_name = "PET_PT", period = "2076_2100")
  
  # Priestley-Taylor potential evapotranspiration with dynamic alpha adjusted
  adjust_PET_annualy(model = M, var_name = "PET_PT", period = "1981_2005", adjustement_period="both")
  adjust_PET_annualy(model = M, var_name = "PET_PT", period = "2076_2100", adjustement_period="both")
  
  # Priestley-Taylor potential evapotranspiration
  assign_ET_var(model = M, var_name = "PET_PT_Zhou", period = "1981_2005")
  assign_ET_var(model = M, var_name = "PET_PT_Zhou", period = "2076_2100")
  
  # Priestley-Taylor potential evapotranspiration adusjted my Morton
  assign_ET_var(model = M, var_name = "PET_PT_Morton", period = "1981_2005")
  assign_ET_var(model = M, var_name = "PET_PT_Morton", period = "2076_2100")
  
  # Makkink_Hansen potential evaporation
  assign_ET_var(model = M, var_name = "PET_Makkink_Hansen", period = "1981_2005")
  assign_ET_var(model = M, var_name = "PET_Makkink_Hansen", period = "2076_2100")
  
  # de Bruin et al. (2016) potential evaporation
  assign_ET_var(model = M, var_name = "PET_deBruin_2016", period = "1981_2005")
  assign_ET_var(model = M, var_name = "PET_deBruin_2016", period = "2076_2100")
  
  # de Bruin and Keijman. (1979) potential evaporation
  assign_ET_var(model = M, var_name = "PET_deBruin_1979", period = "1981_2005")
  assign_ET_var(model = M, var_name = "PET_deBruin_1979", period = "2076_2100")
  
}

################################################################################

##################################
# The deviations from the baseline

# Define the function
get_monthly_perturbation <- function(model, var_name){
  x_baseline <- get(paste0(var_name, "_", model, "_1981_2005"))
  x_future <- get(paste0(var_name, "_", model, "_2076_2100"))
  x <- x_baseline
  x$PERIOD <- "2076_2100-1981_2005"
  
  data_range <- get(paste0(model, "_1981_2005"))
  x[data_range] <- x_future[data_range] - x_baseline[data_range]
  assign(x = paste0("d_", var_name, "_", model) , value = x, envir = parent.frame())
}

# Loop over the models
for(M in models){
  # LE
  get_monthly_perturbation(model = M, var_name = "LE")

  # H
  get_monthly_perturbation(model = M, var_name = "H")
  
  # Rg
  get_monthly_perturbation(model = M, var_name = "Rg")
  
  # SW_net
  get_monthly_perturbation(model = M, var_name = "SW_net")
  
  # LW_net
  get_monthly_perturbation(model = M, var_name = "LW_net")
  
  # Rn
  get_monthly_perturbation(model = M, var_name = "Rn")
  
  # G
  get_monthly_perturbation(model = M, var_name = "G")
  
  # P
  get_monthly_perturbation(model = M, var_name = "P")
  
  # Ta
  get_monthly_perturbation(model = M, var_name = "Ta")
  
  # Saturated vapor pressure
  get_monthly_perturbation(model = M, var_name = "e_sat")
  
  # Actual evapotranspiration
  get_monthly_perturbation(model = M, var_name = "ET")
  
  # Priestley-Taylor potential evapotranspiration with alpha set to 1.26
  get_monthly_perturbation(model = M, var_name = "PET_PT_1.26")
}

################################################################################

######################################################
# The annual (normalized) deviations from the baseline

# Define the function
get_annual_perturbation <- function(model, var_name, fun = "annMeans", normalize = TRUE){
  data_range <- get(paste0(model, "_1981_2005"))[-1] # -1 ensure that the first element of the vector (i.e. ERA5) is skipped
  
  x_baseline <- get(paste0(var_name, "_", model, "_1981_2005"))[,data_range]
  x_delta <- get(paste0("d_", var_name, "_", model))[,data_range]
  
  if(normalize == TRUE){
    if(fun == "annSums"){
      x <- annSums(x_delta)/annSums(x_baseline)
    }else if(fun == "annMeans"){
      x <- annMeans(x_delta)/annMeans(x_baseline)
    }
  }else if(normalize == FALSE){
    if(fun == "annSums"){
      x <- annSums(x_delta)
    }else if(fun == "annMeans"){
      x <- annMeans(x_delta)
    }
  }
  if(normalize == TRUE){
    assign(x = paste0("d_", var_name, "_over_", var_name, "_", model, "_a") , value = x, envir = parent.frame())
  }else if(normalize == FALSE){
    assign(x = paste0("d_", var_name, "_", model, "_a") , value = x, envir = parent.frame())
  }
}

# Loop over the models
for(M in models){
  
  # Precipitation
  get_annual_perturbation(model = M, var_name = "P", fun = "annSums", normalize = TRUE)
  
  # Actual evapotranspiration
  get_annual_perturbation(model = M, var_name = "ET", fun = "annSums", normalize = TRUE)
  
  # Air temperature
  get_annual_perturbation(model = M, var_name = "Ta", fun = "annMeans", normalize = FALSE)
  
  # Saturated vapor pressure
  get_annual_perturbation(model = M, var_name = "e_sat", fun = "annMeans", normalize = TRUE)
  
  # Net radiation
  get_annual_perturbation(model = M, var_name = "Rn", fun = "annMeans", normalize = TRUE)
  
  # Soil heat flux
  get_annual_perturbation(model = M, var_name = "G", fun = "annMeans", normalize = TRUE)
  
  # Global radiation
  get_annual_perturbation(model = M, var_name = "Rg", fun = "annMeans", normalize = TRUE)
  
  # Net shortwave radiation
  get_annual_perturbation(model = M, var_name = "SW_net", fun = "annMeans", normalize = TRUE)
  
  # Net long-wave radiation
  get_annual_perturbation(model = M, var_name = "LW_net", fun = "annMeans", normalize = TRUE)
}

################################################################################

#######################
# The annual statistics

# Define the function
get_annual_stat <- function(model, var_name, period, fun = "annMeans"){
  data_range <- get(paste0(model, "_", period))[-1] # -1 ensure that the first element of the vector (i.e. ERA5) is skipped
  
  x <- get(paste0(var_name, "_", model, "_", period))[,data_range]
 
  if(fun == "annSums"){
    x <- annSums(x)
  }else if(fun == "annMeans"){
    x <- annMeans(x)
  }
  
  assign(x = paste0(var_name, "_", model, "_", period, "_a") , value = x, envir = parent.frame())
}

# Loop over the models
for(M in models){

  # Air temperature
  get_annual_stat(model = M, var_name = "Ta", period = "1981_2005", fun = "annMeans")
  
  # Sensible heat flux
  get_annual_stat(model = M, var_name = "H", period = "1981_2005", fun = "annMeans")
  get_annual_stat(model = M, var_name = "H", period = "2076_2100", fun = "annMeans")
  
  # Latent heat flux
  get_annual_stat(model = M, var_name = "LE", period = "1981_2005", fun = "annMeans")
  get_annual_stat(model = M, var_name = "LE", period = "2076_2100", fun = "annMeans")
  
  # Priestley-Taylor coefficient
  get_annual_stat(model = M, var_name = "alpha_PT", period = "1981_2005", fun = "no")
  
  # Potential evaporation with alpha set to 1.26
  get_annual_stat(model = M, var_name = "PET_PT_1.26", period = "1981_2005", fun = "annSums")
  get_annual_stat(model = M, var_name = "PET_PT_1.26", period = "2076_2100", fun = "annSums")
  
  # Potential evaporation with a model specific (dynamic) alpha
  get_annual_stat(model = M, var_name = "PET_PT", period = "1981_2005", fun = "annSums")
  get_annual_stat(model = M, var_name = "PET_PT", period = "2076_2100", fun = "annSums")
  
  # Potential evaporation following Zhou et al. (2023)
  get_annual_stat(model = M, var_name = "PET_PT_Zhou", period = "1981_2005", fun = "annSums")
  get_annual_stat(model = M, var_name = "PET_PT_Zhou", period = "2076_2100", fun = "annSums")
  
  # Potential evaporation following Morton (1975 and 1976)
  get_annual_stat(model = M, var_name = "PET_PT_Morton", period = "1981_2005", fun = "annSums")
  get_annual_stat(model = M, var_name = "PET_PT_Morton", period = "2076_2100", fun = "annSums")
  
  # Potential evaporation following Hansen (1984)
  get_annual_stat(model = M, var_name = "PET_Makkink_Hansen", period = "1981_2005", fun = "annSums")
  get_annual_stat(model = M, var_name = "PET_Makkink_Hansen", period = "2076_2100", fun = "annSums")
  
  # Potential evaporation following de Bruin (2016)
  get_annual_stat(model = M, var_name = "PET_deBruin_2016", period = "1981_2005", fun = "annSums")
  get_annual_stat(model = M, var_name = "PET_deBruin_2016", period = "2076_2100", fun = "annSums")
  
  # Potential evaporation following de Bruin and Keijman (1979)
  get_annual_stat(model = M, var_name = "PET_deBruin_1979", period = "1981_2005", fun = "annSums")
  get_annual_stat(model = M, var_name = "PET_deBruin_1979", period = "2076_2100", fun = "annSums")
  
  # Precipitation
  get_annual_stat(model = M, var_name = "P", period = "1981_2005", fun = "annSums")
  get_annual_stat(model = M, var_name = "P", period = "2076_2100", fun = "annSums")
  
  # Actual evapotranspiration
  get_annual_stat(model = M, var_name = "ET", period = "1981_2005", fun = "annSums")
  get_annual_stat(model = M, var_name = "ET", period = "2076_2100", fun = "annSums")
}


################################################################################

# Complementary hypothesis, aridity index, evaporative index

# Potential evaporation adjusted by ET based on complementary hypothesis
get_PET_adj_by_CH <- function(model, PET_name, period){

    PET <- get(paste0(PET_name, "_", model, "_", period, "_a"))
    ET <- get(paste0("ET_", model, "_", period, "_a"))
    PET_CH <- 2*PET-ET
    assign(x = paste0(PET_name, "_CH_", model, "_", period, "_a") , value = PET_CH, envir = parent.frame())
}

# Get aridity index and evaporative index
get_Budyko_indicies <- function(model, PET_name, period){
  
  PET <- get(paste0(PET_name, "_", model, "_", period, "_a"))
  P <- get(paste0("P_", model, "_", period, "_a"))
  ET <- get(paste0("ET_", model, "_", period, "_a"))
  
  PET_over_P <- PET/P
  ET_over_P <- ET/P
  
  assign(x = paste0(PET_name, "_over_P_", model, "_", period, "_a") , value = PET_over_P, envir = parent.frame())
  assign(x = paste0("ET_over_P_", model, "_", period, "_a") , value = ET_over_P, envir = parent.frame())
}

# Loop over the models
for(M in models){
  # Potential evaporation adjusted by ET based on complementary hypothesis
  get_PET_adj_by_CH(model = M, PET_name = "PET_PT_1.26", period = "1981_2005")
  get_PET_adj_by_CH(model = M, PET_name = "PET_PT_1.26", period = "2076_2100")
  get_PET_adj_by_CH(model = M, PET_name = "PET_PT", period = "1981_2005")
  get_PET_adj_by_CH(model = M, PET_name = "PET_PT", period = "2076_2100")
  get_PET_adj_by_CH(model = M, PET_name = "PET_PT_Morton", period = "1981_2005")
  get_PET_adj_by_CH(model = M, PET_name = "PET_PT_Morton", period = "2076_2100")
  get_PET_adj_by_CH(model = M, PET_name = "PET_PT_Zhou", period = "1981_2005")
  get_PET_adj_by_CH(model = M, PET_name = "PET_PT_Zhou", period = "2076_2100")
  get_PET_adj_by_CH(model = M, PET_name = "PET_Makkink_Hansen", period = "1981_2005")
  get_PET_adj_by_CH(model = M, PET_name = "PET_Makkink_Hansen", period = "2076_2100")
  get_PET_adj_by_CH(model = M, PET_name = "PET_deBruin_2016", period = "1981_2005")
  get_PET_adj_by_CH(model = M, PET_name = "PET_deBruin_2016", period = "2076_2100")
  get_PET_adj_by_CH(model = M, PET_name = "PET_deBruin_1979", period = "1981_2005")
  get_PET_adj_by_CH(model = M, PET_name = "PET_deBruin_1979", period = "2076_2100")
  
  # Get aridity index and evaporative index
  get_Budyko_indicies(model = M, PET_name = "PET_PT_1.26", period = "1981_2005")
  get_Budyko_indicies(model = M, PET_name = "PET_PT_1.26", period = "2076_2100")
  get_Budyko_indicies(model = M, PET_name = "PET_PT", period = "1981_2005")
  get_Budyko_indicies(model = M, PET_name = "PET_PT", period = "2076_2100")
  get_Budyko_indicies(model = M, PET_name = "PET_PT_Morton", period = "1981_2005")
  get_Budyko_indicies(model = M, PET_name = "PET_PT_Morton", period = "2076_2100")
  get_Budyko_indicies(model = M, PET_name = "PET_PT_Zhou", period = "1981_2005")
  get_Budyko_indicies(model = M, PET_name = "PET_PT_Zhou", period = "2076_2100")
  get_Budyko_indicies(model = M, PET_name = "PET_Makkink_Hansen", period = "1981_2005")
  get_Budyko_indicies(model = M, PET_name = "PET_Makkink_Hansen", period = "2076_2100")
  get_Budyko_indicies(model = M, PET_name = "PET_deBruin_2016", period = "1981_2005")
  get_Budyko_indicies(model = M, PET_name = "PET_deBruin_2016", period = "2076_2100")
  get_Budyko_indicies(model = M, PET_name = "PET_deBruin_1979", period = "1981_2005")
  get_Budyko_indicies(model = M, PET_name = "PET_deBruin_1979", period = "2076_2100")
  
  # Complementary hypothesis
  get_Budyko_indicies(model = M, PET_name = "PET_PT_1.26_CH", period = "1981_2005")
  get_Budyko_indicies(model = M, PET_name = "PET_PT_1.26_CH", period = "2076_2100")
  get_Budyko_indicies(model = M, PET_name = "PET_PT_CH", period = "1981_2005")
  get_Budyko_indicies(model = M, PET_name = "PET_PT_CH", period = "2076_2100")
  get_Budyko_indicies(model = M, PET_name = "PET_PT_Morton_CH", period = "1981_2005")
  get_Budyko_indicies(model = M, PET_name = "PET_PT_Morton_CH", period = "2076_2100")
  get_Budyko_indicies(model = M, PET_name = "PET_PT_Zhou_CH", period = "1981_2005")
  get_Budyko_indicies(model = M, PET_name = "PET_PT_Zhou_CH", period = "2076_2100")
  get_Budyko_indicies(model = M, PET_name = "PET_Makkink_Hansen_CH", period = "1981_2005")
  get_Budyko_indicies(model = M, PET_name = "PET_Makkink_Hansen_CH", period = "2076_2100")
  get_Budyko_indicies(model = M, PET_name = "PET_deBruin_2016_CH", period = "1981_2005")
  get_Budyko_indicies(model = M, PET_name = "PET_deBruin_2016_CH", period = "2076_2100")
  get_Budyko_indicies(model = M, PET_name = "PET_deBruin_1979_CH", period = "1981_2005")
  get_Budyko_indicies(model = M, PET_name = "PET_deBruin_1979_CH", period = "2076_2100")
}
