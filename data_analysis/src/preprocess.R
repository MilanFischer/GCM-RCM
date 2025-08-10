# Unit conversions
W_to_mm <- 0.408*3600*24/10^6

# Psychrometric constant
gamma <- 0.066

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
alpha_PT <- 1.26

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
  y <- 12* colMeans(x,na.rm=TRUE)/mean(c(31,28.25,31,30,31,30,31,31,30,31,30,31))
  
  return(y)
}

# Function for annual sum
annMean <- function(x){
  
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
  y <- colMeans(x,na.rm=TRUE)/mean(c(31,28.25,31,30,31,30,31,31,30,31,30,31))
  
  return(y)
}

######################
# Temporal sub-setting
GCMs_1981_2005_data <- GCMs_data[which(GCMs_data$PERIOD=='1981_2005'),]
GCMs_2076_2100_data <- GCMs_data[which(GCMs_data$PERIOD=='2076_2100'),]

RCMs_1981_2005_data <- RCMs_data[which(RCMs_data$PERIOD=='1981_2005'),]
RCMs_2076_2100_data <- RCMs_data[which(RCMs_data$PERIOD=='2076_2100'),]

# Select the column with variables
GCMs_1981_2005 <- 2:(ncol(GCMs_1981_2005_data)-2)
GCMs_2076_2100 <- 2:(ncol(GCMs_2076_2100_data)-2)
RCMs_1981_2005 <- 2:(ncol(RCMs_1981_2005_data)-2)
RCMs_2076_2100 <- 2:(ncol(RCMs_2076_2100_data)-2)

############################
# Fluxes and state variables

#-----
# GCMs

LE_GCMs_1981_2005 <- GCMs_1981_2005_data[which(GCMs_1981_2005_data$VAR=='LE'),]
LE_GCMs_2076_2100 <- GCMs_2076_2100_data[which(GCMs_2076_2100_data$VAR=='LE'),]

H_GCMs_1981_2005 <- GCMs_1981_2005_data[which(GCMs_1981_2005_data$VAR=='H'),]
H_GCMs_2076_2100 <- GCMs_2076_2100_data[which(GCMs_2076_2100_data$VAR=='H'),]

Rg_GCMs_1981_2005 <- GCMs_1981_2005_data[which(GCMs_1981_2005_data$VAR=='Rg'),]
Rg_GCMs_2076_2100 <- GCMs_2076_2100_data[which(GCMs_2076_2100_data$VAR=='Rg'),]

SW_net_GCMs_1981_2005 <- GCMs_1981_2005_data[which(GCMs_1981_2005_data$VAR=='SW_net'),]
SW_net_GCMs_2076_2100 <- GCMs_2076_2100_data[which(GCMs_2076_2100_data$VAR=='SW_net'),]

LW_net_GCMs_1981_2005 <- GCMs_1981_2005_data[which(GCMs_1981_2005_data$VAR=='LW_net'),]
LW_net_GCMs_2076_2100 <- GCMs_2076_2100_data[which(GCMs_2076_2100_data$VAR=='LW_net'),]

# Correct the sign
LW_net_GCMs_1981_2005[,GCMs_1981_2005] <- (-1)*LW_net_GCMs_1981_2005[,GCMs_1981_2005]
LW_net_GCMs_2076_2100[,GCMs_2076_2100] <- (-1)*LW_net_GCMs_2076_2100[,GCMs_2076_2100]

Rn_GCMs_1981_2005 <- SW_net_GCMs_1981_2005
Rn_GCMs_1981_2005[,GCMs_1981_2005] <- SW_net_GCMs_1981_2005[,GCMs_1981_2005]+LW_net_GCMs_1981_2005[,GCMs_1981_2005]
Rn_GCMs_2076_2100 <- SW_net_GCMs_2076_2100
Rn_GCMs_2076_2100[,GCMs_1981_2005] <- SW_net_GCMs_2076_2100[,GCMs_2076_2100]+LW_net_GCMs_2076_2100[,GCMs_2076_2100]
Rn_GCMs_1981_2005$VAR <- 'Rn'
Rn_GCMs_2076_2100$VAR <- 'Rn'

G_GCMs_1981_2005 <- Rn_GCMs_1981_2005
G_GCMs_1981_2005[,GCMs_1981_2005] <- Rn_GCMs_1981_2005[,GCMs_1981_2005]-H_GCMs_1981_2005[,GCMs_1981_2005]-LE_GCMs_1981_2005[,GCMs_1981_2005]
G_GCMs_2076_2100 <- Rn_GCMs_2076_2100
G_GCMs_2076_2100[,GCMs_2076_2100] <- Rn_GCMs_2076_2100[,GCMs_2076_2100]-H_GCMs_2076_2100[,GCMs_2076_2100]-LE_GCMs_2076_2100[,GCMs_2076_2100]
G_GCMs_1981_2005$VAR <- 'G'
G_GCMs_2076_2100$VAR <- 'G'

P_GCMs_1981_2005 <- GCMs_1981_2005_data[which(GCMs_1981_2005_data$VAR=='P'),]
P_GCMs_2076_2100 <- GCMs_2076_2100_data[which(GCMs_2076_2100_data$VAR=='P'),]

Ta_GCMs_1981_2005 <- GCMs_1981_2005_data[which(GCMs_1981_2005_data$VAR=='Ta'),]
Ta_GCMs_2076_2100 <- GCMs_2076_2100_data[which(GCMs_2076_2100_data$VAR=='Ta'),]

n_days_GCMs <-matrix(rep(c(31,28.25,31,30,31,30,31,31,30,31,30,31),ncol(GCMs_1981_2005_data)-4),nrow=12,ncol=ncol(GCMs_1981_2005_data)-4)
ET_GCMs_1981_2005 <- LE_GCMs_1981_2005
ET_GCMs_1981_2005[,GCMs_1981_2005] <- LE_GCMs_1981_2005[,GCMs_1981_2005]*n_days_GCMs*W_to_mm
ET_GCMs_2076_2100 <- LE_GCMs_2076_2100
ET_GCMs_2076_2100[,GCMs_2076_2100] <- LE_GCMs_2076_2100[,GCMs_2076_2100]*n_days_GCMs*W_to_mm
ET_GCMs_2076_2100 <- LE_GCMs_2076_2100


# Clausius-Clapeyron
e_sat_GCMs_1981_2005 <- Ta_GCMs_1981_2005
e_sat_GCMs_2076_2100 <- Ta_GCMs_2076_2100

e_sat_GCMs_1981_2005[,GCMs_1981_2005] <- a*exp(b*Ta_GCMs_1981_2005[,GCMs_1981_2005]/(Ta_GCMs_1981_2005[,GCMs_1981_2005]+c))
e_sat_GCMs_2076_2100[,GCMs_2076_2100] <- a*exp(b*Ta_GCMs_2076_2100[,GCMs_2076_2100]/(Ta_GCMs_2076_2100[,GCMs_2076_2100]+c))
e_sat_GCMs_1981_2005$VAR <- 'e_sat'
e_sat_GCMs_2076_2100$VAR <- 'e_sat'


SVP_GCMs_1981_2005 <- Ta_GCMs_1981_2005
SVP_GCMs_2076_2100 <- Ta_GCMs_2076_2100

SVP_GCMs_1981_2005[,GCMs_1981_2005] <- b*c* e_sat_GCMs_1981_2005[,GCMs_1981_2005]/(c + Ta_GCMs_1981_2005[,GCMs_1981_2005])^2
SVP_GCMs_2076_2100[,GCMs_2076_2100] <- b*c* e_sat_GCMs_2076_2100[,GCMs_2076_2100]/(c + Ta_GCMs_2076_2100[,GCMs_2076_2100])^2
SVP_GCMs_1981_2005$VAR <- 'SVP'
SVP_GCMs_2076_2100$VAR <- 'SVP'

# Equilibrium evaporation
ET_eq_GCMs_1981_2005 <- Ta_GCMs_1981_2005
ET_eq_GCMs_2076_2100 <- Ta_GCMs_2076_2100

ET_eq_GCMs_1981_2005[,GCMs_1981_2005] <- SVP_GCMs_1981_2005[,GCMs_1981_2005]*(Rn_GCMs_1981_2005[,GCMs_1981_2005]-G_GCMs_1981_2005[,GCMs_1981_2005])/(SVP_GCMs_1981_2005[,GCMs_1981_2005] + gamma)*n_days_GCMs*W_to_mm
ET_eq_GCMs_2076_2100[,GCMs_2076_2100] <- SVP_GCMs_2076_2100[,GCMs_2076_2100]*(Rn_GCMs_2076_2100[,GCMs_2076_2100]-G_GCMs_2076_2100[,GCMs_2076_2100])/(SVP_GCMs_2076_2100[,GCMs_2076_2100] + gamma)*n_days_GCMs*W_to_mm

ET_eq_GCMs_1981_2005$VAR <- 'Equilibrium_ET'
ET_eq_GCMs_2076_2100$VAR <- 'Equilibrium_ET'

# Replace negative values by zero
ET_eq_GCMs_1981_2005[,GCMs_1981_2005][ET_eq_GCMs_1981_2005[,GCMs_1981_2005]<0] <- 0
ET_eq_GCMs_2076_2100[,GCMs_2076_2100][ET_eq_GCMs_2076_2100[,GCMs_2076_2100]<0] <- 0

# Priestley-Taylor potential evapotranspiration with alpha set to 1.26
PET_PT_1.26_GCMs_1981_2005 <- Ta_GCMs_1981_2005
PET_PT_1.26_GCMs_2076_2100 <- Ta_GCMs_2076_2100

PET_PT_1.26_GCMs_1981_2005[,GCMs_1981_2005] <- alpha_PT*ET_eq_GCMs_1981_2005[,GCMs_1981_2005]
PET_PT_1.26_GCMs_2076_2100[,GCMs_2076_2100] <- alpha_PT*ET_eq_GCMs_2076_2100[,GCMs_2076_2100]

PET_PT_1.26_GCMs_1981_2005$VAR <- 'PET_PT_alpha_1.26'
PET_PT_1.26_GCMs_2076_2100$VAR <- 'PET_PT_alpha_1.26'

# Dynamic alpha
alpha_PT_GCMs_1981_2005 <- Ta_GCMs_1981_2005

alpha_PT_GCMs_1981_2005[,GCMs_1981_2005] <- ET_GCMs_1981_2005[,GCMs_1981_2005]/ET_eq_GCMs_1981_2005[,GCMs_1981_2005]

alpha_PT_GCMs_1981_2005 <- alpha_PT_GCMs_1981_2005[alpha_PT_GCMs_1981_2005$MONTH %in% c('IV','V','VI','VII','VIII'),]
tmp <- apply(alpha_PT_GCMs_1981_2005[,GCMs_1981_2005],2,max)

alpha_PT_GCMs_1981_2005 <- alpha_PT_GCMs_1981_2005[1,]
alpha_PT_GCMs_1981_2005[,GCMs_1981_2005] <- tmp

alpha_PT_GCMs_1981_2005$MONTH <- 'IV-VIII'
alpha_PT_GCMs_1981_2005$VAR <- 'alpha_PT'

PET_PT_GCMs_1981_2005 <- PET_PT_1.26_GCMs_1981_2005
PET_PT_GCMs_1981_2005$VAR <- 'PET_PT_alpha_dynamic'

for(i in GCMs_1981_2005){
  PET_PT_GCMs_1981_2005[,i] <- ET_eq_GCMs_1981_2005[,i]*max(alpha_PT_GCMs_1981_2005[i],alpha_PT)
}

# Check that annual sum of ET is not exceeding PET and correct it if yes
ET_GCMs_a_1981_2005 <- annSums(ET_GCMs_1981_2005[,GCMs_1981_2005])
PET_GCMs_a_1981_2005 <- annSums(PET_PT_GCMs_1981_2005[,GCMs_1981_2005])
ET_eq_GCMs_a_1981_2005 <- annSums(ET_eq_GCMs_1981_2005[,GCMs_1981_2005])

tmp <- alpha_PT_GCMs_1981_2005
tmp[,GCMs_1981_2005] <- ET_GCMs_a_1981_2005/ET_eq_GCMs_a_1981_2005

for(i in GCMs_1981_2005){
  alpha_PT_GCMs_1981_2005[,i] <- ifelse(tmp[,i]>alpha_PT&tmp[,i]>alpha_PT_GCMs_1981_2005[,i],
                                        tmp[,i], alpha_PT_GCMs_1981_2005[,i])
  PET_PT_GCMs_1981_2005[,i] <- ET_eq_GCMs_1981_2005[,i]*max(alpha_PT_GCMs_1981_2005[i],alpha_PT)
}
rm(list=c('ET_GCMs_a_1981_2005','PET_GCMs_a_1981_2005','ET_eq_GCMs_a_1981_2005','tmp'))

PET_PT_GCMs_2076_2100 <- PET_PT_1.26_GCMs_2076_2100
PET_PT_GCMs_2076_2100$VAR <- 'PET_PT_alpha_dynamic'

for(i in GCMs_2076_2100){
  PET_PT_GCMs_2076_2100[,i] <- ET_eq_GCMs_2076_2100[,i]*max(alpha_PT_GCMs_1981_2005[i],alpha_PT)
}

#-----
# RCMs

LE_RCMs_1981_2005 <- RCMs_1981_2005_data[which(RCMs_1981_2005_data$VAR=='LE'),]
LE_RCMs_2076_2100 <- RCMs_2076_2100_data[which(RCMs_2076_2100_data$VAR=='LE'),]

H_RCMs_1981_2005 <- RCMs_1981_2005_data[which(RCMs_1981_2005_data$VAR=='H'),]
H_RCMs_2076_2100 <- RCMs_2076_2100_data[which(RCMs_2076_2100_data$VAR=='H'),]

Rg_RCMs_1981_2005 <- RCMs_1981_2005_data[which(RCMs_1981_2005_data$VAR=='Rg'),]
Rg_RCMs_2076_2100 <- RCMs_2076_2100_data[which(RCMs_2076_2100_data$VAR=='Rg'),]

SW_net_RCMs_1981_2005 <- RCMs_1981_2005_data[which(RCMs_1981_2005_data$VAR=='SW_net'),]
SW_net_RCMs_2076_2100 <- RCMs_2076_2100_data[which(RCMs_2076_2100_data$VAR=='SW_net'),]

LW_net_RCMs_1981_2005 <- RCMs_1981_2005_data[which(RCMs_1981_2005_data$VAR=='LW_net'),]
LW_net_RCMs_2076_2100 <- RCMs_2076_2100_data[which(RCMs_2076_2100_data$VAR=='LW_net'),]

# Correct the sign
LW_net_RCMs_1981_2005[,RCMs_1981_2005] <- (-1)*LW_net_RCMs_1981_2005[,RCMs_1981_2005]
LW_net_RCMs_2076_2100[,RCMs_2076_2100] <- (-1)*LW_net_RCMs_2076_2100[,RCMs_2076_2100]

Rn_RCMs_1981_2005 <- SW_net_RCMs_1981_2005
Rn_RCMs_1981_2005[,RCMs_1981_2005] <- SW_net_RCMs_1981_2005[,RCMs_1981_2005]+LW_net_RCMs_1981_2005[,RCMs_1981_2005]
Rn_RCMs_2076_2100 <- SW_net_RCMs_2076_2100
Rn_RCMs_2076_2100[,RCMs_1981_2005] <- SW_net_RCMs_2076_2100[,RCMs_2076_2100]+LW_net_RCMs_2076_2100[,RCMs_2076_2100]
Rn_RCMs_1981_2005$VAR <- 'Rn'
Rn_RCMs_2076_2100$VAR <- 'Rn'

G_RCMs_1981_2005 <- Rn_RCMs_1981_2005
G_RCMs_1981_2005[,RCMs_1981_2005] <- Rn_RCMs_1981_2005[,RCMs_1981_2005]-H_RCMs_1981_2005[,RCMs_1981_2005]-LE_RCMs_1981_2005[,RCMs_1981_2005]
G_RCMs_2076_2100 <- Rn_RCMs_2076_2100
G_RCMs_2076_2100[,RCMs_2076_2100] <- Rn_RCMs_2076_2100[,RCMs_2076_2100]-H_RCMs_2076_2100[,RCMs_2076_2100]-LE_RCMs_2076_2100[,RCMs_2076_2100]
G_RCMs_1981_2005$VAR <- 'G'
G_RCMs_2076_2100$VAR <- 'G'

P_RCMs_1981_2005 <- RCMs_1981_2005_data[which(RCMs_1981_2005_data$VAR=='P'),]
P_RCMs_2076_2100 <- RCMs_2076_2100_data[which(RCMs_2076_2100_data$VAR=='P'),]

Ta_RCMs_1981_2005 <- RCMs_1981_2005_data[which(RCMs_1981_2005_data$VAR=='Ta'),]
Ta_RCMs_2076_2100 <- RCMs_2076_2100_data[which(RCMs_2076_2100_data$VAR=='Ta'),]

n_days_RCMs <- matrix(rep(c(31,28.25,31,30,31,30,31,31,30,31,30,31),ncol(RCMs_1981_2005_data)-4),nrow=12,ncol=ncol(RCMs_1981_2005_data)-4)
ET_RCMs_1981_2005 <- LE_RCMs_1981_2005
ET_RCMs_1981_2005[,RCMs_1981_2005] <- LE_RCMs_1981_2005[,RCMs_1981_2005]*n_days_RCMs*W_to_mm
ET_RCMs_2076_2100 <- LE_RCMs_2076_2100
ET_RCMs_2076_2100[,RCMs_2076_2100] <- LE_RCMs_2076_2100[,RCMs_2076_2100]*n_days_RCMs*W_to_mm
ET_RCMs_2076_2100 <- LE_RCMs_2076_2100

# Clausius-Clapeyron
e_sat_RCMs_1981_2005 <- Ta_RCMs_1981_2005
e_sat_RCMs_2076_2100 <- Ta_RCMs_2076_2100

e_sat_RCMs_1981_2005[,RCMs_1981_2005] <- a*exp(b*Ta_RCMs_1981_2005[,RCMs_1981_2005]/(Ta_RCMs_1981_2005[,RCMs_1981_2005]+c))
e_sat_RCMs_2076_2100[,RCMs_2076_2100] <- a*exp(b*Ta_RCMs_2076_2100[,RCMs_2076_2100]/(Ta_RCMs_2076_2100[,RCMs_2076_2100]+c))

e_sat_RCMs_1981_2005$VAR <- 'e_sat'
e_sat_RCMs_2076_2100$VAR <- 'e_sat'

SVP_RCMs_1981_2005 <- Ta_RCMs_1981_2005
SVP_RCMs_2076_2100 <- Ta_RCMs_2076_2100

SVP_RCMs_1981_2005[,RCMs_1981_2005] <- b*c* e_sat_RCMs_1981_2005[,RCMs_1981_2005]/(c + Ta_RCMs_1981_2005[,RCMs_1981_2005])^2
SVP_RCMs_2076_2100[,RCMs_2076_2100] <- b*c* e_sat_RCMs_2076_2100[,RCMs_2076_2100]/(c + Ta_RCMs_2076_2100[,RCMs_2076_2100])^2

SVP_RCMs_1981_2005$VAR <- 'SVP'
SVP_RCMs_2076_2100$VAR <- 'SVP'

# Equilibrium evaporation
ET_eq_RCMs_1981_2005 <- Ta_RCMs_1981_2005
ET_eq_RCMs_2076_2100 <- Ta_RCMs_2076_2100

ET_eq_RCMs_1981_2005[,RCMs_1981_2005] <- SVP_RCMs_1981_2005[,RCMs_1981_2005]*(Rn_RCMs_1981_2005[,RCMs_1981_2005]-G_RCMs_1981_2005[,RCMs_1981_2005])/(SVP_RCMs_1981_2005[,RCMs_1981_2005] + gamma)*n_days_RCMs*W_to_mm
ET_eq_RCMs_2076_2100[,RCMs_2076_2100] <- SVP_RCMs_2076_2100[,RCMs_2076_2100]*(Rn_RCMs_2076_2100[,RCMs_2076_2100]-G_RCMs_2076_2100[,RCMs_2076_2100])/(SVP_RCMs_2076_2100[,RCMs_2076_2100] + gamma)*n_days_RCMs*W_to_mm

ET_eq_RCMs_1981_2005$VAR <- 'Equilibrium_ET'
ET_eq_RCMs_2076_2100$VAR <- 'Equilibrium_ET'

# Replace negative values by zero
ET_eq_RCMs_1981_2005[,RCMs_1981_2005][ET_eq_RCMs_1981_2005[,RCMs_1981_2005]<0] <- 0
ET_eq_RCMs_2076_2100[,RCMs_2076_2100][ET_eq_RCMs_2076_2100[,RCMs_2076_2100]<0] <- 0

# Priestley-Taylor potential evapotranspiration with alpha set to 1.26
PET_PT_1.26_RCMs_1981_2005 <- Ta_RCMs_1981_2005
PET_PT_1.26_RCMs_2076_2100 <- Ta_RCMs_2076_2100

PET_PT_1.26_RCMs_1981_2005[,RCMs_1981_2005] <- alpha_PT*ET_eq_RCMs_1981_2005[,RCMs_1981_2005]
PET_PT_1.26_RCMs_2076_2100[,RCMs_2076_2100] <- alpha_PT*ET_eq_RCMs_2076_2100[,RCMs_2076_2100]

PET_PT_1.26_RCMs_1981_2005$VAR <- 'PET_PT_alpha_1.26'
PET_PT_1.26_RCMs_2076_2100$VAR <- 'PET_PT_alpha_1.26'

# Dynamic alpha
alpha_PT_RCMs_1981_2005 <- Ta_RCMs_1981_2005

alpha_PT_RCMs_1981_2005[,RCMs_1981_2005] <- ET_RCMs_1981_2005[,RCMs_1981_2005]/ET_eq_RCMs_1981_2005[,RCMs_1981_2005]

alpha_PT_RCMs_1981_2005 <- alpha_PT_RCMs_1981_2005[alpha_PT_RCMs_1981_2005$MONTH %in% c('IV','V','VI','VII','VIII'),]
tmp <- apply(alpha_PT_RCMs_1981_2005[,RCMs_1981_2005],2,max)

alpha_PT_RCMs_1981_2005 <- alpha_PT_RCMs_1981_2005[1,]
alpha_PT_RCMs_1981_2005[,RCMs_1981_2005] <- tmp

alpha_PT_RCMs_1981_2005$MONTH <- 'IV-VIII'
alpha_PT_RCMs_1981_2005$VAR <- 'alpha_PT'

PET_PT_RCMs_1981_2005 <- PET_PT_1.26_RCMs_1981_2005
PET_PT_RCMs_1981_2005$VAR <- 'PET_PT_alpha_dynamic'

for(i in RCMs_1981_2005){
  PET_PT_RCMs_1981_2005[,i] <- ET_eq_RCMs_1981_2005[,i]*max(alpha_PT_RCMs_1981_2005[i],alpha_PT)
}

# Check that annual sum of ET is not exceeding PET and correct it if yes
ET_RCMs_a_1981_2005 <- annSums(ET_RCMs_1981_2005[,RCMs_1981_2005])
PET_RCMs_a_1981_2005 <- annSums(PET_PT_RCMs_1981_2005[,RCMs_1981_2005])
ET_eq_RCMs_a_1981_2005 <- annSums(ET_eq_RCMs_1981_2005[,RCMs_1981_2005])

tmp <- alpha_PT_RCMs_1981_2005
tmp[,RCMs_1981_2005] <- ET_RCMs_a_1981_2005/ET_eq_RCMs_a_1981_2005

for(i in RCMs_1981_2005){
  alpha_PT_RCMs_1981_2005[,i] <- ifelse(tmp[,i]>alpha_PT&tmp[,i]>alpha_PT_RCMs_1981_2005[,i],
                                        tmp[,i], alpha_PT_RCMs_1981_2005[,i])
  PET_PT_RCMs_1981_2005[,i] <- ET_eq_RCMs_1981_2005[,i]*max(alpha_PT_RCMs_1981_2005[i],alpha_PT)
}
rm(list=c('ET_RCMs_a_1981_2005','PET_RCMs_a_1981_2005','ET_eq_RCMs_a_1981_2005','tmp'))

PET_PT_RCMs_2076_2100 <- PET_PT_1.26_RCMs_2076_2100
PET_PT_RCMs_2076_2100$VAR <- 'PET_PT_alpha_dynamic'

for(i in RCMs_2076_2100){
  PET_PT_RCMs_2076_2100[,i] <- ET_eq_RCMs_2076_2100[,i]*max(alpha_PT_RCMs_1981_2005[i],alpha_PT)
}

##################################
# The deviations from the baseline

#-----
# GCMs

d_LE_GCMs <- LE_GCMs_2076_2100; d_LE_GCMs$PERIOD <- '2076_2100-1981_2005'
d_LE_GCMs[,GCMs_2076_2100] <- LE_GCMs_2076_2100[,GCMs_2076_2100]-LE_GCMs_1981_2005[,GCMs_2076_2100]
d_H_GCMs <- H_GCMs_2076_2100; H_GCMs_2076_2100$PERIOD <- '2076_2100-1981_2005'
d_H_GCMs[,GCMs_2076_2100] <- d_H_GCMs[,GCMs_2076_2100]-H_GCMs_1981_2005[,GCMs_2076_2100]

d_Rg_GCMs <- Rg_GCMs_2076_2100; d_Rg_GCMs$PERIOD <- '2076_2100-1981_2005'
d_Rg_GCMs[,GCMs_2076_2100] <- Rg_GCMs_2076_2100[,GCMs_2076_2100]-Rg_GCMs_1981_2005[,GCMs_2076_2100]
d_SW_net_GCMs <- SW_net_GCMs_2076_2100; d_SW_net_GCMs$PERIOD <- '2076_2100-1981_2005'
d_SW_net_GCMs[,GCMs_2076_2100] <- SW_net_GCMs_2076_2100[,GCMs_2076_2100]-SW_net_GCMs_1981_2005[,GCMs_2076_2100]
d_LW_net_GCMs <- LW_net_GCMs_2076_2100; d_LW_net_GCMs$PERIOD <- '2076_2100-1981_2005'
d_LW_net_GCMs[,GCMs_2076_2100] <- LW_net_GCMs_2076_2100[,GCMs_2076_2100]-LW_net_GCMs_1981_2005[,GCMs_2076_2100]
d_Rn_GCMs <- Rn_GCMs_2076_2100; d_Rn_GCMs$PERIOD <- '2076_2100-1981_2005'
d_Rn_GCMs[,GCMs_2076_2100] <- Rn_GCMs_2076_2100[,GCMs_2076_2100]-Rn_GCMs_1981_2005[,GCMs_2076_2100]

d_P_GCMs <- P_GCMs_2076_2100; d_P_GCMs$PERIOD <- '2076_2100-1981_2005'
d_P_GCMs[,GCMs_2076_2100] <- P_GCMs_2076_2100[,GCMs_2076_2100]-P_GCMs_1981_2005[,GCMs_2076_2100]
d_Ta_GCMs <- Ta_GCMs_2076_2100; d_Ta_GCMs$PERIOD <- '2076_2100-1981_2005'
d_Ta_GCMs[,GCMs_2076_2100] <- Ta_GCMs_2076_2100[,GCMs_2076_2100]-Ta_GCMs_1981_2005[,GCMs_2076_2100]

d_e_sat_GCMs <- e_sat_GCMs_2076_2100; d_e_sat_GCMs$PERIOD <- '2076_2100-1981_2005'
d_e_sat_GCMs[,GCMs_2076_2100] <- e_sat_GCMs_2076_2100[,GCMs_2076_2100]-e_sat_GCMs_1981_2005[,GCMs_2076_2100]

d_ET_GCMs <- ET_GCMs_2076_2100; d_ET_GCMs$PERIOD <- '2076_2100-1981_2005'
d_ET_GCMs[,GCMs_2076_2100] <- ET_GCMs_2076_2100[,GCMs_2076_2100]-ET_GCMs_1981_2005[,GCMs_2076_2100]

d_PET_PT_1.26_GCMs <- PET_PT_1.26_GCMs_2076_2100; d_PET_PT_1.26_GCMs$PERIOD <- '2076_2100-1981_2005'
d_PET_PT_1.26_GCMs[,GCMs_2076_2100] <- PET_PT_1.26_GCMs_2076_2100[,GCMs_2076_2100]-PET_PT_1.26_GCMs_1981_2005[,GCMs_2076_2100]

#-----
# RCMs

d_LE_RCMs <- LE_RCMs_2076_2100; d_LE_RCMs$PERIOD <- '2076_2100-1981_2005'
d_LE_RCMs[,RCMs_2076_2100] <- LE_RCMs_2076_2100[,RCMs_2076_2100]-LE_RCMs_1981_2005[,RCMs_2076_2100]
d_H_RCMs <- H_RCMs_2076_2100; H_RCMs_2076_2100$PERIOD <- '2076_2100-1981_2005'
d_H_RCMs[,RCMs_2076_2100] <- d_H_RCMs[,RCMs_2076_2100]-H_RCMs_1981_2005[,RCMs_2076_2100]

d_Rg_RCMs <- Rg_RCMs_2076_2100; d_Rg_RCMs$PERIOD <- '2076_2100-1981_2005'
d_Rg_RCMs[,RCMs_2076_2100] <- Rg_RCMs_2076_2100[,RCMs_2076_2100]-Rg_RCMs_1981_2005[,RCMs_2076_2100]
d_SW_net_RCMs <- SW_net_RCMs_2076_2100; d_SW_net_RCMs$PERIOD <- '2076_2100-1981_2005'
d_SW_net_RCMs[,RCMs_2076_2100] <- SW_net_RCMs_2076_2100[,RCMs_2076_2100]-SW_net_RCMs_1981_2005[,RCMs_2076_2100]
d_LW_net_RCMs <- LW_net_RCMs_2076_2100; d_LW_net_RCMs$PERIOD <- '2076_2100-1981_2005'
d_LW_net_RCMs[,RCMs_2076_2100] <- LW_net_RCMs_2076_2100[,RCMs_2076_2100]-LW_net_RCMs_1981_2005[,RCMs_2076_2100]
d_Rn_RCMs <- Rn_RCMs_2076_2100; d_Rn_RCMs$PERIOD <- '2076_2100-1981_2005'
d_Rn_RCMs[,RCMs_2076_2100] <- Rn_RCMs_2076_2100[,RCMs_2076_2100]-Rn_RCMs_1981_2005[,RCMs_2076_2100]

d_P_RCMs <- P_RCMs_2076_2100; d_P_RCMs$PERIOD <- '2076_2100-1981_2005'
d_P_RCMs[,RCMs_2076_2100] <- P_RCMs_2076_2100[,RCMs_2076_2100]-P_RCMs_1981_2005[,RCMs_2076_2100]
d_Ta_RCMs <- Ta_RCMs_2076_2100; d_Ta_RCMs$PERIOD <- '2076_2100-1981_2005'
d_Ta_RCMs[,RCMs_2076_2100] <- Ta_RCMs_2076_2100[,RCMs_2076_2100]-Ta_RCMs_1981_2005[,RCMs_2076_2100]

d_e_sat_RCMs <- e_sat_RCMs_2076_2100; d_e_sat_RCMs$PERIOD <- '2076_2100-1981_2005'
d_e_sat_RCMs[,RCMs_2076_2100] <- e_sat_RCMs_2076_2100[,RCMs_2076_2100]-e_sat_RCMs_1981_2005[,RCMs_2076_2100]

d_ET_RCMs <- ET_RCMs_2076_2100; d_ET_RCMs$PERIOD <- '2076_2100-1981_2005'
d_ET_RCMs[,RCMs_2076_2100] <- ET_RCMs_2076_2100[,RCMs_2076_2100]-ET_RCMs_1981_2005[,RCMs_2076_2100]

d_PET_PT_1.26_RCMs <- PET_PT_1.26_RCMs_2076_2100; d_PET_PT_1.26_RCMs$PERIOD <- '2076_2100-1981_2005'
d_PET_PT_1.26_RCMs[,RCMs_2076_2100] <- PET_PT_1.26_RCMs_2076_2100[,RCMs_2076_2100]-PET_PT_1.26_RCMs_1981_2005[,RCMs_2076_2100]

