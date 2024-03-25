# An R script to compute pFDR and the  variance 
# -> Integral_pFDR (alpha, beta, m0, m1) # computes pFDR
# -> pFDR_variance (alpha, beta, m0, m1) # computes the  variance 

#  	    -- where alpha is Probability of the type I error, 
#	    		 (1-beta) is Probability of type II error,
#	    		 m0 is the number of truly null hypotheses, 
#	    		 and m1 is the number of alternative hypotheses

#Note
#===================================================
#The False Discovery Rate and its Variance Integral 
#Formulas have been implemented in R and thoroughly tested with varying values of alpha, beta, m0, and m1. 
#In most numerical simulation scenarios, the functions return exact values. However, it has been observed that 
#for large values of alpha, m0, and m1, the code may return negative values due to insufficient grid density 
#for integration. This issue can be resolved by incorporating a finer integral grid


#calculates the pFDR using the integration approach
pFDR_IntegralApproach <- function(alpha, beta, m0, m1) {
  #Calculates the pFDR function at a given value during integration
  integrant <- function(s, alpha, beta, m0, m1) {
    result = (s * alpha  + 1 - alpha)^(m0 - 1) * (s * beta + 1 - beta)^m1
    return(result)
  }
  
  #these are the ranges where the integrations are explicited calculated in
  bounds <- c(0, 0.0001, 0.0005, 0.001, 0.01, 0.025, 0.05, 
              0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 
              0.95, 0.96, 0.97, 0.98, 0.985, 0.99, 0.995, 1)
  pFDR_theoretical = 0 #value of the pFDR_theoretical
  
  for(i in 2:length(bounds)) { #calculate each range of integration sequential and combine
    lower = bounds[i-1] #lower bound of integration
    upper = bounds[i] #upper bound of integration
    pFDR_theoretical = pFDR_theoretical + integrate(integrant, lower=lower, upper=upper, 
                                                    alpha=alpha, beta=beta, m0=m0, m1=m1)$value
  }
  
  normalization <- (m0 * alpha) / (1 - (1-alpha)^m0 * (1-beta)^m1)
  pFDR_theoretical <-  pFDR_theoretical * normalization #apply normalization value   
  return (pFDR_theoretical);
}
