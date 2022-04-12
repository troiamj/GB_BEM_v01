###############################################################################################################################
###############################################################################################################################
size_at_age <- 6.382417                                           # from Groeschel-Taylor et al
resamp001_W2 <- size_at_age                                       # grams wet mass
prey_ED <- 3698.0                                                 # joules per gram of wet mass
oxycal_coeff <- 13560.0                                           # joules per gram of oxygen

###############################################################################################################################
###############################################################################################################################
#CONSUMPTION function                  # equation 2
consumption2 <- function(T,               # environmental temperature (degrees C) at which consumption will be calculated
                         p,               # proportion of maximum consumption
                         CA_eq2,          # intercept of allometric mass function
                         CB_eq2,          # slope of allometric mass function; should be negative bc bigger fish consume less per gram of body mass than smaller fish
                         W,               # weight of fish (grams)
                         CTM_eq2,         # critical thermal maximum (degrees C)
                         CTO_eq2,         # laboratory temperature preferendum (degrees C)
                         CQ_eq2)          # approximates a Q10; the rate at which the function increases over relatively low water temperatures
{Y <- log(CQ_eq2)*(CTM_eq2-CTO_eq2+2)
Z <- log(CQ_eq2)*(CTM_eq2-CTO_eq2)
X <- (Z^2*(1+(1+40/Y)^0.5)^2)/400
V <- (CTM_eq2-T)/(CTM_eq2-CTO_eq2)
fT_C <- V^X*exp(X*(1-V))
Cmax <- CA_eq2*W^CB_eq2
C <- Cmax*p*fT_C                         
return(C)}                             # function returns the specific consumption rate (grams of food consumed per gram of fish mass per day)

###############################################################################################################################
###############################################################################################################################
# RESPIRATION function                 # equation 1
respiration1 <- function(T,
                         ACT,
                         RA,
                         RB,
                         W,
                         RQ)
{
  fT <- exp(RQ*T)
  ACTIVITY <- ACT
  R <- RA*W^RB*fT*ACTIVITY
  return(R)
}

###############################################################################################################################
###############################################################################################################################
#EGESTION function
# no function for this species

###############################################################################################################################
###############################################################################################################################
#EXCRETION function
# no function for this species

###############################################################################################################################
###############################################################################################################################
# read in BEM parameters
BEM_parameters <- read.csv("<<YOUR DIRECTORY HERE>>\\troia_etal_BEM_v01\\input_bem_parameters_v01\\BEM_parameters_v01.csv", header = TRUE)
# parameter set to use for the model run --> resamp001
BEM_parameters <- BEM_parameters[which(BEM_parameters$genspe=="resamp001"),]

###############################################################################################################################
###############################################################################################################################
# make parms dataframe to store simulated bioenergetics parameters for each day of a 365-day year

# read in modeled temperature
# put each column (COMID) into its own list element
modeled_WT_HUC12100201 <- read.csv("<<YOUR DIRECTORY HERE>>\\troia_etal_BEM_v01\\input_temperature_v01\\modeled_WT_HUC12100201.csv", header = TRUE)
# reformat so day 1 (row 1) is not julian day 1 (jan 1st), but rather julian day 60 (first day of spring months)
#modeled_WT_HUC12100201 <- rbind(modeled_WT_HUC12100201[121:365,], modeled_WT_HUC12100201[1:120,])

# rename headers; replace "X" with "COMID"
new_header <- gsub("X", "COMID_", names(modeled_WT_HUC12100201))
colnames(modeled_WT_HUC12100201) <- new_header

# put each column (COMID) into its own list element
list.sim_parms <- list()
for(k in 1:ncol(modeled_WT_HUC12100201)){list.sim_parms[[k]] <- modeled_WT_HUC12100201[,k]}
names(list.sim_parms) <- colnames(modeled_WT_HUC12100201)

# add julian date
for(k in 1:ncol(modeled_WT_HUC12100201)){list.sim_parms[[k]] <- data.frame(c(1:365), list.sim_parms[[k]])}
for(k in 1:ncol(modeled_WT_HUC12100201)){colnames(list.sim_parms[[k]]) <- c("julian","temp")}
# add parm columns with zeros as fillers
for(k in 1:ncol(modeled_WT_HUC12100201)){list.sim_parms[[k]]$C1_ins <- NA}
for(k in 1:ncol(modeled_WT_HUC12100201)){list.sim_parms[[k]]$C2_ins <- NA}
for(k in 1:ncol(modeled_WT_HUC12100201)){list.sim_parms[[k]]$R1_ins <- NA}
for(k in 1:ncol(modeled_WT_HUC12100201)){list.sim_parms[[k]]$R2_ins <- NA}
for(k in 1:ncol(modeled_WT_HUC12100201)){list.sim_parms[[k]]$F1_ins <- NA}
for(k in 1:ncol(modeled_WT_HUC12100201)){list.sim_parms[[k]]$F2_ins <- NA}
for(k in 1:ncol(modeled_WT_HUC12100201)){list.sim_parms[[k]]$U1_ins <- NA}
for(k in 1:ncol(modeled_WT_HUC12100201)){list.sim_parms[[k]]$U2_ins <- NA}
for(k in 1:ncol(modeled_WT_HUC12100201)){list.sim_parms[[k]]$SDA1_ins <- NA}
for(k in 1:ncol(modeled_WT_HUC12100201)){list.sim_parms[[k]]$SDA2_ins <- NA}
for(k in 1:ncol(modeled_WT_HUC12100201)){list.sim_parms[[k]]$W1_ins <- NA}
for(k in 1:ncol(modeled_WT_HUC12100201)){list.sim_parms[[k]]$W2_ins <- NA}
for(k in 1:ncol(modeled_WT_HUC12100201)){list.sim_parms[[k]]$C1_cum <- NA}
for(k in 1:ncol(modeled_WT_HUC12100201)){list.sim_parms[[k]]$C2_cum <- NA}
for(k in 1:ncol(modeled_WT_HUC12100201)){list.sim_parms[[k]]$R1_cum <- NA}
for(k in 1:ncol(modeled_WT_HUC12100201)){list.sim_parms[[k]]$R2_cum <- NA}
for(k in 1:ncol(modeled_WT_HUC12100201)){list.sim_parms[[k]]$F1_cum <- NA}
for(k in 1:ncol(modeled_WT_HUC12100201)){list.sim_parms[[k]]$F2_cum <- NA}
for(k in 1:ncol(modeled_WT_HUC12100201)){list.sim_parms[[k]]$U1_cum <- NA}
for(k in 1:ncol(modeled_WT_HUC12100201)){list.sim_parms[[k]]$U2_cum <- NA}
for(k in 1:ncol(modeled_WT_HUC12100201)){list.sim_parms[[k]]$SDA1_cum <- NA}
for(k in 1:ncol(modeled_WT_HUC12100201)){list.sim_parms[[k]]$SDA2_cum <- NA}
for(k in 1:ncol(modeled_WT_HUC12100201)){list.sim_parms[[k]]$W1_cum <- NA}
for(k in 1:ncol(modeled_WT_HUC12100201)){list.sim_parms[[k]]$W2_cum <- NA}
# add starting weight in first row; add zeros to all other first rows
for(k in 1:ncol(modeled_WT_HUC12100201)){list.sim_parms[[k]][1,c("C1_ins")] <- 0}
for(k in 1:ncol(modeled_WT_HUC12100201)){list.sim_parms[[k]][1,c("C2_ins")] <- 0}
for(k in 1:ncol(modeled_WT_HUC12100201)){list.sim_parms[[k]][1,c("R1_ins")] <- 0}
for(k in 1:ncol(modeled_WT_HUC12100201)){list.sim_parms[[k]][1,c("R2_ins")] <- 0}
for(k in 1:ncol(modeled_WT_HUC12100201)){list.sim_parms[[k]][1,c("F1_ins")] <- 0}
for(k in 1:ncol(modeled_WT_HUC12100201)){list.sim_parms[[k]][1,c("F2_ins")] <- 0}
for(k in 1:ncol(modeled_WT_HUC12100201)){list.sim_parms[[k]][1,c("U1_ins")] <- 0}
for(k in 1:ncol(modeled_WT_HUC12100201)){list.sim_parms[[k]][1,c("U2_ins")] <- 0}
for(k in 1:ncol(modeled_WT_HUC12100201)){list.sim_parms[[k]][1,c("SDA1_ins")] <- 0}
for(k in 1:ncol(modeled_WT_HUC12100201)){list.sim_parms[[k]][1,c("SDA2_ins")] <- 0}
for(k in 1:ncol(modeled_WT_HUC12100201)){list.sim_parms[[k]][1,c("W1_ins")] <- 0}
for(k in 1:ncol(modeled_WT_HUC12100201)){list.sim_parms[[k]][1,c("W2_ins")] <- 0}
for(k in 1:ncol(modeled_WT_HUC12100201)){list.sim_parms[[k]][1,c("C1_cum")] <- 0}
for(k in 1:ncol(modeled_WT_HUC12100201)){list.sim_parms[[k]][1,c("C2_cum")] <- 0}
for(k in 1:ncol(modeled_WT_HUC12100201)){list.sim_parms[[k]][1,c("R1_cum")] <- 0}
for(k in 1:ncol(modeled_WT_HUC12100201)){list.sim_parms[[k]][1,c("R2_cum")] <- 0}
for(k in 1:ncol(modeled_WT_HUC12100201)){list.sim_parms[[k]][1,c("F1_cum")] <- 0}
for(k in 1:ncol(modeled_WT_HUC12100201)){list.sim_parms[[k]][1,c("F2_cum")] <- 0}
for(k in 1:ncol(modeled_WT_HUC12100201)){list.sim_parms[[k]][1,c("U1_cum")] <- 0}
for(k in 1:ncol(modeled_WT_HUC12100201)){list.sim_parms[[k]][1,c("U2_cum")] <- 0}
for(k in 1:ncol(modeled_WT_HUC12100201)){list.sim_parms[[k]][1,c("SDA1_cum")] <- 0}
for(k in 1:ncol(modeled_WT_HUC12100201)){list.sim_parms[[k]][1,c("SDA2_cum")] <- 0}
for(k in 1:ncol(modeled_WT_HUC12100201)){list.sim_parms[[k]][1,c("W1_cum")] <- 0}
for(k in 1:ncol(modeled_WT_HUC12100201)){list.sim_parms[[k]][1,c("W2_cum")] <- resamp001_W2}  # set initial weight to age1

# sort columns alphabetically
modeled_WT_HUC12100201 <- modeled_WT_HUC12100201[,order(colnames(modeled_WT_HUC12100201))]

###############################################################################################################################
###############################################################################################################################
for(j in 1:length(list.sim_parms)){
  for(i in 2:nrow(list.sim_parms[[j]])){
    # simulate consumption --> grams of prey
    list.sim_parms[[j]][i,c("C1_ins")] <- ifelse(is.nan(consumption2(list.sim_parms[[j]][i,c("temp")],     # temperature on the ith julian day
                                                                     BEM_parameters$CP,                                  
                                                                     BEM_parameters$CA,
                                                                     BEM_parameters$CB,
                                                                     list.sim_parms[[j]][i-1,c("W2_cum")], # mass on the preceding day
                                                                     BEM_parameters$CTM,
                                                                     BEM_parameters$CTO,
                                                                     BEM_parameters$CQ)),
                                                 0,                                                        # ** if temp > CTM, then consumption is zero
                                                 consumption2(list.sim_parms[[j]][i,c("temp")],            # temperature on the ith julian day
                                                              BEM_parameters$CP,                                         
                                                              BEM_parameters$CA,
                                                              BEM_parameters$CB,
                                                              list.sim_parms[[j]][i-1,c("W2_cum")],        # mass on the preceding day
                                                              BEM_parameters$CTM,
                                                              BEM_parameters$CTO,
                                                              BEM_parameters$CQ)) * list.sim_parms[[j]][i-1,c("W2_cum")]
    
    # simulate consumption --> joules of energy
    list.sim_parms[[j]][i,c("C2_ins")] <- list.sim_parms[[j]][i,c("C1_ins")] * prey_ED                     # convert from grams of food to joules with prey energy density parameter
    
    # simulate respiration --> grams of oxygen
    list.sim_parms[[j]][i,c("R1_ins")] <- ifelse(is.nan(respiration1(list.sim_parms[[j]][i,c("temp")],     # temperature on the ith julian day
                                                                     BEM_parameters$ACT,                                  
                                                                     BEM_parameters$RA,
                                                                     BEM_parameters$RB,
                                                                     list.sim_parms[[j]][i-1,c("W2_cum")],  # mass on the preceding day
                                                                     BEM_parameters$RQ)),
                                                 NA,                                                       # ** if temp > RTM (lethal limit), then the fish theoretically dies, so assign NA value
                                                 respiration1(list.sim_parms[[j]][i,c("temp")],            # temperature on the ith julian day
                                                              BEM_parameters$ACT,                                         
                                                              BEM_parameters$RA,
                                                              BEM_parameters$RB,
                                                              list.sim_parms[[j]][i-1,c("W2_cum")],         # mass on the preceding day
                                                              BEM_parameters$RQ)) * list.sim_parms[[j]][i-1,c("W2_cum")]
    
    # simulate respiration --> joules of energy
    list.sim_parms[[j]][i,c("R2_ins")] <- list.sim_parms[[j]][i,c("R1_ins")] * oxycal_coeff                      # convert from grams of oxygen to joules with oxycalorific coefficient
    
    # simulate egestion --> grams of prey
    list.sim_parms[[j]][i,c("F1_ins")] <- list.sim_parms[[j]][i,c("C1_ins")] * BEM_parameters$UA                    # assume egestion is a constant proportion of consumption, for now...
    
    # simulate egestion --> joules of energy
    list.sim_parms[[j]][i,c("F2_ins")] <- list.sim_parms[[j]][i,c("C2_ins")] * BEM_parameters$FA                    # assume excretion is a constant proportion of consumption, for now...
    
    # simulate excretion --> grams of prey
    list.sim_parms[[j]][i,c("U1_ins")] <- list.sim_parms[[j]][i,c("C1_ins")] * BEM_parameters$UA                    # assume excretion is a constant proportion of consumption, for now...  
    
    # simulate excretion --> joules of energy
    list.sim_parms[[j]][i,c("U2_ins")] <- list.sim_parms[[j]][i,c("C2_ins")] * BEM_parameters$UA                    # assume excretion is a constant proportion of consumption, for now...  
    
    # simulate specific dynamic action --> grams of prey
    list.sim_parms[[j]][i,c("SDA1_ins")] <- BEM_parameters$SDA * (list.sim_parms[[j]][i,c("C1_ins")] - list.sim_parms[[j]][i,c("F1_ins")])      # assume SDA is a constant proportion of assimilated energy (consumption minus egestion)  
    
    # simulate specific dynamic action --> joules of energy
    list.sim_parms[[j]][i,c("SDA2_ins")] <- BEM_parameters$SDA * (list.sim_parms[[j]][i,c("C2_ins")] - list.sim_parms[[j]][i,c("F2_ins")])      # assume SDA is a constant proportion of assimilated energy (consumption minus egestion)  
    
    # simulate daily weight change, that is: C-(M+E+U+SDA) --> joules of energy
    list.sim_parms[[j]][i,c("W1_ins")] <- list.sim_parms[[j]][i,c("C2_ins")] - (list.sim_parms[[j]][i,c("R2_ins")] + list.sim_parms[[j]][i,c("F2_ins")] + list.sim_parms[[j]][i,c("U2_ins")] + list.sim_parms[[j]][i,c("SDA2_ins")])
    
    # simulate daily weight change, that is: C-(M+E+U+SDA) --> grams of body mass
    list.sim_parms[[j]][i,c("W2_ins")] <- list.sim_parms[[j]][i,c("W1_ins")] * (1/BEM_parameters$ED)                                                         # convert from joules to grams of body mass with predator energy density parameter
    
    # simulate cumulative weight --> joules of energy
    list.sim_parms[[j]][i,c("W1_cum")] <- list.sim_parms[[j]][i-1,c("W1_cum")] + list.sim_parms[[j]][i,c("W1_ins")]
    
    # simulate cumulative weight --> grams of body mass
    list.sim_parms[[j]][i,c("W2_cum")] <- list.sim_parms[[j]][i-1,c("W2_cum")] + list.sim_parms[[j]][i,c("W2_ins")]
  }
  
  # compute remaining cumulative parmeters
  list.sim_parms[[j]]$C1_cum <- cumsum(list.sim_parms[[j]]$C1_ins)
  list.sim_parms[[j]]$C2_cum <- cumsum(list.sim_parms[[j]]$C2_ins)
  list.sim_parms[[j]]$R1_cum <- cumsum(list.sim_parms[[j]]$R1_ins)
  list.sim_parms[[j]]$R2_cum <- cumsum(list.sim_parms[[j]]$R2_ins)
  list.sim_parms[[j]]$F1_cum <- cumsum(list.sim_parms[[j]]$F1_ins)
  list.sim_parms[[j]]$F2_cum <- cumsum(list.sim_parms[[j]]$F2_ins)
  list.sim_parms[[j]]$U1_cum <- cumsum(list.sim_parms[[j]]$U1_ins)
  list.sim_parms[[j]]$U2_cum <- cumsum(list.sim_parms[[j]]$U2_ins)
  list.sim_parms[[j]]$SDA1_cum <- cumsum(list.sim_parms[[j]]$SDA1_ins)
  list.sim_parms[[j]]$SDA2_cum <- cumsum(list.sim_parms[[j]]$SDA2_ins)
}

###############################################################################################################################
# add COMID column
for(i in 1:length(list.sim_parms)){list.sim_parms[[i]]$COMID <- names(list.sim_parms)[i]}

# put all reaches in one dataframe
sim_parms <- do.call("rbind", list.sim_parms)
rownames(sim_parms) <- NULL

###############################################################################################################################
# write to dbf
write.csv(sim_parms, "<<YOUR DIRECTORY HERE>>\\troia_etal_BEM_v01\\output_bem_projections_v01\\bem_HUC12100201.csv", row.names = FALSE)
