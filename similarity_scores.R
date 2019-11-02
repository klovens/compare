#Coefficient of Quartile Deviation
#$$\text{Coefficient of Quartile Deviation = }\frac{Q_3 - Q_1}{Q_3 + Q_1} \times 100$$

calculate_quartile_dev <- function(exp_mat){
  coefficients <- c()
  for(seq in 1:nrow(exp_mat)){
    Q <- quantile(as.numeric(exp_mat[seq,])) 
    coeff_qd <- (Q[4] - Q[2])/(Q[4] + Q[2])
    coefficients <- c(coefficients, coeff_qd*100)
  }
  return(coefficients)
}

coefficient_variation <- function(exp_mat){
  coefficients <- c()
  for(seq in 1:nrow(exp_mat)){
    ex <- as.numeric(exp_mat[seq,])
    mu <- mean(ex)
    s <- sd(ex)
    coefficients <- c(coefficients, (s/mu)*100)
  }
  return(coefficients)
}

