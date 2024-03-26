# An R script to compute pFDR and the  variance 
# -> pFDR_IntegralApproach (alpha, beta, m0, m1) # computes pFDR
# -> FDPvariance_IntegralApproach (alpha, beta, m0, m1) # computes the  variance 

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


#Calculates the function value for integration over two variables for the first term of the FDP Variance equation 
firstTermIntegrant <- function(x1, x2, alpha, beta, m0, m1) {
  temp = (x1 * x2 * alpha  + 1 - alpha)^(m0 - 1) * (x1 * x2 * beta + 1 - beta)^m1
  return (temp);
}

#Calculates the 1d integrated value over a single variable using the oned_integrant function for the first term of the FDP Variance equation 
firstTermIntegral <- function(x2, alpha, beta, m0, m1) {
  boundSet1 <- c(0, 0.7, 0.8, 0.9, #set of ranges to calculate 1d integration of the first term where the function is less volatile
                 0.95, 0.96, 0.97, 0.98, 0.99)
  boundSet2 <- c(0.99, 0.9998, 1) #set of ranges to calculate 1d integration at higher subdivision count where the function is more volatile
  result = 0 #value of 1d integration
  
  for(i in 2:length(boundSet1)) { #sum 1d integration values at lower subdivision count
    lower = boundSet1[i-1]
    upper = boundSet1[i]
    result = result + integrate(firstTermIntegrant, lower=lower, upper=upper, x2=x2, 
                                alpha=alpha, beta=beta, m0=m0, m1=m1)$value
  }
  
  for(i in 2:length(boundSet2)) { #sum 1d integration values at higher subdivision count where the function is more volatile
    lower = boundSet2[i-1]
    upper = boundSet2[i]
    result = result + integrate(firstTermIntegrant, lower=lower, upper=upper, x2=x2, 
                                alpha=alpha, beta=beta, m0=m0, m1=m1, subdivisions=100000)$value
  }
  
  return (result)
}

#calculate the first term of the Variance of the FDP formula 
firstTermValue <- function(alpha, beta, m0, m1) {
  boundSet1 <- c(0, 0.7, 0.8, 0.9, 0.925, 0.95, 0.96, 0.97, 0.98, 0.99) #set of ranges to calculate 2d integration of the first term where the function is less volatile
  boundSet2 <- c(0.99, 0.9998, 1) #set of ranges to calculate 2d integration at higher subdivision count where the function is more volatile
  result = 0 #value of the first term
  
  for(i in 2:length(boundSet1)) { #sum 2d integration values at lower subdivision count
    lower = boundSet1[i-1]
    upper = boundSet1[i]
    result = result + integrate(Vectorize(firstTermIntegral), 
                                alpha=alpha, beta=beta, m0=m0, m1=m1, lower=lower, upper=upper)$value
  }
  
  for(i in 2:length(boundSet2)) { #sum 2d integration values at higher subdivision count where the function is more volatile
    lower = boundSet2[i-1]
    upper = boundSet2[i]
    result = result + integrate(Vectorize(firstTermIntegral), 
                                alpha=alpha, beta=beta, m0=m0, m1=m1, lower=lower, upper=upper, subdivisions = 10000)$value
  }
  
  normalization <- (m0 * alpha) / (1 - (1-alpha)^m0 * (1-beta)^m1)
  normalized = result * normalization #apply normalization constant to first term
  return (normalized) 
}


#Calculates the function value for integration over two variables for the second term of the FDP Variance equation 
secondTermIntegrant <- function(x1, x2, alpha, beta, m0, m1) {
  temp = x1 * x2 * 
    (x1 * x2 * alpha  + 1 - alpha)^(m0 - 2) * 
    (x1 * x2 * beta + 1 - beta)^m1
  return (temp);
}

#Calculates the 1d integrated value over a single variable using the oned_2_integrant function for the second term of the FDP Variance equation 
secondTermIntegral <- function(x2, alpha, beta, m0, m1) {
  boundSet1 <- c(0, 0.6, 0.7, 0.8, 0.9, 0.91, #set of ranges to calculate 1d integration of the second term where the function is less volatile
                 0.95, 0.96, 0.97, 0.98, 0.99, 0.995)
  boundSet2 <- c(0.995, 0.9998, 1) #set of ranges to calculate 1d integration at higher subdivision count where the function is more volatile
  result = 0 #value of 1d integration
  
  for(i in 2:length(boundSet1)) { #sum 1d integration values at lower subdivision count
    lower = boundSet1[i-1]
    upper = boundSet1[i]
    result = result + integrate(secondTermIntegrant, lower=lower, upper=upper, x2=x2, 
                                alpha=alpha, beta=beta, m0=m0, m1=m1)$value
  }
  
  for(i in 2:length(boundSet2)) { #sum 1d integration values at higher subdivision count where the function is more volatile
    lower = boundSet2[i-1]
    upper = boundSet2[i]
    result = result + integrate(secondTermIntegrant, lower=lower, upper=upper, x2=x2, 
                                alpha=alpha, beta=beta, m0=m0, m1=m1, subdivisions = 10000)$value
  }
  
  return (result)
}

#calculate the second term of the Variance of the FDP formula 
secondTermValue <- function(alpha, beta, m0, m1) {
  bounds <- c(0, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.925, 
              0.95, 0.96, 0.97, 0.98, 0.99, 0.9998, 1) #set of ranges to calculate 2d integration of the second term
  result = 0 #value of the second term
  
  for(i in 2:length(bounds)) { #sum 2d integration values
    lower = bounds[i-1]
    upper = bounds[i]
    result = result + integrate(Vectorize(secondTermIntegral), 
                                alpha=alpha, beta=beta, m0=m0, m1=m1, lower=lower, upper=upper)$value
  }
  
  normalization <- (m0 * (m0- 1) * alpha^2) / (1 - (1 - alpha)^m0*(1 - beta)^m1)
  normalized = result * normalization #apply normalization constant to second term
  return (normalized) 
}


#Computes the first term of the variance equation: For large values of m0 & m1 (in millions) and small alpha & beta integration ranges should be smaller 
thirdTermSplitValue <- function(alpha, beta, m0, m1) {
  length = 10000;
  x = seq(0, 1, len = length);
  return_value = 0.000;
  
  for(i in 1:length) {
    temp = (x[i] * alpha + 1 - alpha)^(m0-1) * (x[i] * beta + 1 - beta)^m1;
    
    if(temp > 0.0000001) {
      break;
    }
  }
  
  if(i < length && i > 1) {
    return_value = x[i-1];
  }
  return (return_value);
}

#Calculates the pFDR function at a given value during integration
thirdTermIntegrant <- function(s, alpha, beta, m0, m1)  {
  result = (s * alpha  + 1 - alpha)^(m0 - 1) * (s * beta + 1 - beta)^m1
  return(result)
}

#calculate the third term of the Variance of the FDP formula 
thirdTermValue <- function(alpha, beta, m0, m1) {
  bounds <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9,
              0.95, 0.96, 0.97, 0.98, 0.99) #set of ranges to calculate 1d integration of the third term at lower subdivision count
  result = 0 #value of the third term
  split_value = thirdTermSplitValue(alpha, beta, m0, m1); #determine half way point, and calculate third term in two portins
  
  # The script first determines the first value of the argument (x) for which the integrand is larger than 10^(-7). 
  # Then it does the integration for two ranges: up to the determined value and after that.
  if(split_value > 0.0) { 
    result <- integrate(thirdTermIntegrant, lower=0, upper=split_value, 
                        alpha=alpha, beta=beta, m0=m0, m1=m1)$value +
      integrate(thirdTermIntegrant, lower=split_value, upper=1, 
                alpha=alpha, beta=beta, m0=m0, m1=m1)$value
    
  } else {
    for(i in 2:length(bounds)) { #sum 1d integration values at lower subdivision count
      lower = bounds[i-1]
      upper = bounds[i]
      result = result + integrate(thirdTermIntegrant, lower=lower, upper=upper, 
                                  alpha=alpha, beta=beta, m0=m0, m1=m1)$value
    }
    
    result = (result +
                integrate(thirdTermIntegrant, lower=0.99, upper=0.9998, 
                          alpha=alpha, beta=beta, m0=m0, m1=m1, subdivisions=10000)$value +
                integrate(thirdTermIntegrant, lower=0.9998, upper=1, 
                          alpha=alpha, beta=beta, m0=m0, m1=m1)$value
    )  #sum 1d integration values at higher subdivision count where the function is more volatile
  }
  
  normalization = (m0 * alpha) / (1 - (1-alpha)^m0 * (1-beta)^m1)
  normalized = (result * normalization)^2 #for the variance, the square is needed
  return (normalized) 
}


#calculates the variance of the FDP
FDPvariance_IntegralApproach <- function(alpha, beta, m0, m1) {# the variance is the made of the three terms
  
  term1 = firstTermValue(alpha, beta, m0, m1) #compute first term
  term2 = secondTermValue(alpha, beta, m0, m1) #compute second term
  term3 = thirdTermValue(alpha, beta, m0, m1) #compute third term
  var = term1 + term2 - term3 #combine terms
  
  if(var < 0.0) { #if variance is negative throw warning
    print("The integral calculations went wrong: Need to increase the grid density");
    return (-1);
  }
  
  options(digits= 5);
  return(var); #return calculated variance
}