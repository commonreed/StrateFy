# LA: leaf area mm^2
# LFW: leaf fresh weight mg
# LDW: leaf dry weight mg

# LA
# LDMC (leaf dry matter content) = LDW/LFW * 100
# SLA (specific leaf area) = LA/LDW

strateFy <- function(Name, LA, LDMC, SLA){
  # transformations
  LA_1 <- sqrt(LA/894205) * 100
  LDMC_1 <- log((LDMC/100)/(1-(LDMC/100)))
  SLA_1 <- log(SLA)
  
  # Regression against calibration PCA
  PC2_C <- -0.8678 + 1.6464 * LA_1
  PC1_S <- 1.3369+0.000010019*(1-exp(-0.0000000000022303*LDMC_1))+4.5835*(1-exp(-0.2328*LDMC_1))
  PC1_R <- -57.5924 + 62.6802*exp(-0.0288*SLA_1)
  
  # CSR classification boundary definition
  C_min <- 0
  S_min <- -0.756451214853076
  R_min <- -11.3467682227961
  
  #Negetive outlier correction CSR
  C_neg <- PC2_C; C_neg[C_neg < C_min] <- C_min; C_neg
  S_neg <- PC1_S; S_neg[S_neg < S_min] <- S_min; S_neg
  R_neg <- PC1_R; R_neg[R_neg < R_min] <- R_min; R_neg
  
  C_max <- 57.3756711966087
  S_max <- 5.79158377609218
  R_max <-1.10795515716546
  
  # Positive outlier correction Csr
  C_pos <- C_neg; C_pos[C_pos > C_max] <- C_max; C_pos
  S_pos <- S_neg; S_pos[S_pos > S_max] <- S_max; S_pos
  R_pos <- R_neg; R_pos[R_pos > R_max] <- R_max; R_pos
  
  # +ve translation (values ranges spanning zero are shifted to all become positive)
  # Positive translation coef. CSR
  C_coef <- abs(C_min)
  S_coef <- abs(S_min)
  R_coef <- abs(R_min)
  
  C_pos_t <- C_pos + C_coef
  S_pos_t <- S_pos + S_coef
  R_pos_t <- R_pos + R_coef
  
  # Range PCA positive translation
  PC2_C_range <- C_max + abs(C_min)
  PC1_S_range <- S_max + abs(S_min)
  PC1_R_range <- R_max + abs(R_min)
  
  # Proportion of total variability (CSR)
  C_pro <- C_pos_t/PC2_C_range * 100
  S_pro <- S_pos_t/PC1_S_range * 100
  R_pro <- 100 - R_pos_t/PC1_R_range * 100
  
  # Percentage Conversion Coefficient
  Per_coeff <- 100/(C_pro + S_pro + R_pro)
  
  C <- C_pro * Per_coeff
  S <- S_pro * Per_coeff
  R <- R_pro * Per_coeff
  
  res <- data.frame(Name, C, S, R)
  
  
  CSR_class <- c('C', 'C/CR', 'C/CS', 'CR', 'C/CSR',' CS', 'CR/CSR', 'CS/CSR', 
                 'R/CR', 'CSR', 'S/CS', 'R/CSR', 'S/CSR', 'R', 'SR/CSR', 
                 'S',' R/SR', 'S/SR',' SR')
  C_threshold <- c(90.000, 72.500, 72.500, 47.500, 54.200, 47.500, 41.675, 
                   41.675, 22.500, 100/3, 22.500, 22.900, 22.900, 5.000, 
                   16.650, 5.000, 5.000, 5.000, 5.000)
  S_threshold <- c(5.000, 5.000, 22.500, 5.000, 22.900, 47.500, 16.650, 
                   41.675, 5.000, 100/3, 72.500, 22.900, 54.200, 5.000, 
                   41.675, 90.000, 22.500, 72.500, 47.500)
  R_threshold <- c(5.000, 22.500, 5.000, 47.500, 22.900, 5.000, 41.675,
                   16.650, 72.500, 100/3, 5.000, 54.200, 22.900, 90.000,
                   41.675, 5.000, 72.500, 22.500, 47.500)
  
  for(i in 1:length(Name)){
    variance <- (C_threshold - C[i])^2 + (S_threshold - S[i])^2 + (R_threshold - R[i])
    res[i, "CSR_Type"] <- CSR_class[which.min(variance)]
  }
  return(res)
}


#Name <- c("a", "b", "c")
#LA <- c(11, 55, 36061)
#LDMC <- c(29, 13, 35)
#SLA <- c(6, 34, 14)

#calCSR(Name, LA, LDMC, SLA)
