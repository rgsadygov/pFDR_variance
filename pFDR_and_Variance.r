
#
#
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
#Formulas have been implemented in R and thoroughly 
#tested with varying values of alpha, beta, m0, and m1. 
#In most numerical simulation scenarios, the functions 
#return exact values. However, it has been observed that 
#for large values of alpha, m0, and m1, the code may 
#return negative values due to insufficient grid density 
#for integration. This issue can be resolved by incorporating 
#a finer integral grid, as demonstrated on lines 40 -62, 127 - 147, 213 - 237, 309 - 324




Integral_pFDR <- function(alpha, beta, m0, m1)
{   
  integrant <- function(s, alpha = 0.1, beta = 0.1, m0 = 1, m1 = 1) 
  {
    
    (s * alpha  + 1 - alpha)^(m0 - 1) *
      (s * beta + 1 - beta)^m1 ;
  }
  
  pFDR_theoretical <- integrate(integrant, lower=0, upper = 0.0001, alpha = alpha, beta = beta, m0 = m0, m1 = m1)$value +
    
    integrate(integrant, lower=0.0001, upper = .0005, alpha = alpha, beta = beta, m0 = m0, m1 = m1)$value +
    integrate(integrant, lower=0.0005, upper = .001, alpha = alpha, beta = beta, m0 = m0, m1 = m1)$value +
    integrate(integrant, lower=0.001, upper = .01, alpha = alpha, beta = beta, m0 = m0, m1 = m1)$value +
    integrate(integrant, lower=0.01, upper = .025, alpha = alpha, beta = beta, m0 = m0, m1 = m1)$value +
    integrate(integrant, lower=0.025, upper = .05, alpha = alpha, beta = beta, m0 = m0, m1 = m1)$value +
    integrate(integrant, lower=0.05, upper = .1, alpha = alpha, beta = beta, m0 = m0, m1 = m1)$value +
    integrate(integrant, lower=0.1, upper = .2, alpha = alpha, beta = beta, m0 = m0, m1 = m1)$value +
    integrate(integrant, lower=0.2, upper = .30, alpha = alpha, beta = beta, m0 = m0, m1 = m1)$value +
    integrate(integrant, lower=0.3, upper = .40, alpha = alpha, beta = beta, m0 = m0, m1 = m1)$value +
    integrate(integrant, lower=0.4, upper = .50, alpha = alpha, beta = beta, m0 = m0, m1 = m1)$value +
    integrate(integrant, lower=0.5, upper = .60, alpha = alpha, beta = beta, m0 = m0, m1 = m1)$value +
    
    integrate(integrant, lower=0.6, upper = .70, alpha = alpha, beta = beta, m0 = m0, m1 = m1)$value +
    integrate(integrant, lower=0.70, upper = .80, alpha = alpha, beta = beta, m0 = m0, m1 = m1)$value +
    integrate(integrant, lower=0.80, upper = .90, alpha = alpha, beta = beta, m0 = m0, m1 = m1)$value +
    integrate(integrant, lower=0.9, upper = .95, alpha = alpha, beta = beta, m0 = m0, m1 = m1)$value +
    integrate(integrant, lower=0.95, upper = 0.96, alpha = alpha, beta = beta, m0 = m0, m1 = m1)$value +
    integrate(integrant, lower=0.96, upper = 0.97, alpha = alpha, beta = beta, m0 = m0, m1 = m1)$value +
    integrate(integrant, lower=0.97, upper = 0.98, alpha = alpha, beta = beta, m0 = m0, m1 = m1)$value +
    integrate(integrant, lower=0.98, upper = 0.985, alpha = alpha, beta = beta, m0 = m0, m1 = m1)$value +
    integrate(integrant, lower=0.985, upper = 0.99, alpha = alpha, beta = beta, m0 = m0, m1 = m1)$value +
    integrate(integrant, lower=0.99, upper = 0.995, alpha = alpha, beta = beta, m0 = m0, m1 = m1)$value +
    integrate(integrant, lower=0.995, upper = 1., alpha = alpha, beta = beta, m0 = m0, m1 = m1)$value
  
  #print(pFDR_theoretical);
  FDR <- pFDR_theoretical * m0 * alpha;
  pFDR_theoretical <-  pFDR_theoretical * m0 * alpha / (1. - (1. - alpha)^m0 * (1. - beta)^m1 );  
  return (pFDR_theoretical);
}


return_limit <- function(alpha, beta, m0, m1)
{
#
# this function computes the first term of the three terms
# in the expression for the variance:
#
# For large values of m0 and m1 (in millions) and small alpha and bete
# the ranges of integrations should de made shorter. 
#
#
  length = 10000;
  
  x = seq(0, 1, len = length);
  
  return_value = 0.000;
  
  for(i in 1:length)
  {
    temp = (x[i] * alpha + 1 - alpha)^(m0-1) * (x[i] * beta + 1 - beta)^m1;
    
    if(temp > 0.0000001)
    {
      break;
    }
    
  }
  
  if(i < length && i > 1)
  {
    return_value = x[i-1];
  }
  
  #  print(return_value);
  #  print((x[1] * alpha + 1 - alpha)^(m0-1) * (x[1] * beta + 1 - beta)^m1 );
  
  return (return_value);
}



first_term <- function(alpha, beta, m0, m1)
{
  
  oned_integrant <- function(x1, x2, alpha, beta, m0, m1) 
  {
    
    temp = (x1 * x2 * alpha  + 1 - alpha)^(m0 - 1) *
      (x1 * x2 * beta + 1 - beta)^m1
    
    return (temp);
  }
  
  
  
  oned_integral <- function(x2, alpha, beta, m0, m1)
  {
    res = integrate(oned_integrant, lower=0, upper= 0.7, x2 = x2, alpha=alpha, beta = beta, m0 = m0, m1 = m1)$value
    
    res = res + integrate(oned_integrant, lower=0.7, upper= 0.8, x2 = x2, alpha=alpha, beta = beta, m0 = m0, m1 = m1)$value
    
    res = res + integrate(oned_integrant, lower=0.8, upper= 0.9, x2 = x2, alpha=alpha, beta = beta, m0 = m0, m1 = m1)$value
    
    res = res + integrate(oned_integrant, lower=0.9, upper= 0.95, x2 = x2, alpha=alpha, beta = beta, m0 = m0, m1 = m1)$value
    
    res = res + integrate(oned_integrant, lower=0.95, upper= 0.96, x2 = x2, alpha=alpha, beta = beta, m0 = m0, m1 = m1)$value
    
    res = res + integrate(oned_integrant, lower=0.96, upper= 0.97, x2 = x2, alpha=alpha, beta = beta, m0 = m0, m1 = m1)$value
    
    res = res + integrate(oned_integrant, lower=0.97, upper= 0.98, x2 = x2, alpha=alpha, beta = beta, m0 = m0, m1 = m1)$value
    
    res = res + integrate(oned_integrant, lower=0.98, upper= 0.99, x2 = x2, alpha=alpha, beta = beta, m0 = m0, m1 = m1)$value
    
    #res = res + integrate(oned_integrant, lower=0.99, upper= 1.0, x2 = x2, alpha=alpha, beta = beta, m0 = m0, m1 = m1)$value
    
    res = res + integrate(oned_integrant, lower=0.99, upper= 0.9998, x2 = x2, alpha=alpha, beta = beta, m0 = m0, m1 = m1, subdivisions = 100000)$value
    
    res = res + integrate(oned_integrant, lower=0.9998, upper= 1.0, x2 = x2, alpha=alpha, beta = beta, m0 = m0, m1 = m1, subdivisions = 100000)$value
    
    return (res)
  }
  
  
  
  twod = integrate(Vectorize(oned_integral), alpha=alpha, beta=beta, m0=m0, m1=m1, lower=0, upper=0.7)$value
  
  twod = twod + integrate(Vectorize(oned_integral), alpha=alpha, beta=beta, m0=m0, m1=m1, lower=0.7, upper=0.8)$value
  
  twod = twod + integrate(Vectorize(oned_integral), alpha=alpha, beta=beta, m0=m0, m1=m1, lower=0.8, upper=0.9)$value
  
  twod = twod + integrate(Vectorize(oned_integral), alpha=alpha, beta=beta, m0=m0, m1=m1, lower=0.9, upper=0.925)$value
  
  twod = twod + integrate(Vectorize(oned_integral), alpha=alpha, beta=beta, m0=m0, m1=m1, lower=0.925, upper=0.95)$value
  
  twod = twod + integrate(Vectorize(oned_integral), alpha=alpha, beta=beta, m0=m0, m1=m1, lower=0.95, upper=0.96)$value
  
  twod = twod + integrate(Vectorize(oned_integral), alpha=alpha, beta=beta, m0=m0, m1=m1, lower=0.96, upper=0.97)$value
  
  twod = twod + integrate(Vectorize(oned_integral), alpha=alpha, beta=beta, m0=m0, m1=m1, lower=0.97, upper=0.98)$value
  
  twod = twod + integrate(Vectorize(oned_integral), alpha=alpha, beta=beta, m0=m0, m1=m1, lower=0.98, upper=0.99)$value
  
  #twod = twod + integrate(Vectorize(oned_integral), alpha=alpha, beta=beta, m0=m0, m1=m1, lower=0.99, upper=1.)$value
  
  twod = twod + integrate(Vectorize(oned_integral), alpha=alpha, beta=beta, m0=m0, m1=m1, lower=0.99, upper=0.9998, subdivisions=10000)$value
  
  twod = twod + integrate(Vectorize(oned_integral), alpha=alpha, beta=beta, m0=m0, m1=m1, lower=0.9998, upper=1., subdivisions=10000)$value
  
  
  
  # print(c(integrate(Vectorize(oned_integral), alpha=alpha, beta=beta, m0=m0, m1=m1, lower=0., upper=1., subdivisions=1000000)$value,
  #     twod));
  
  twod = twod * m0 * alpha  / (1 - (1 - alpha)^m0*(1 - beta)^m1);
  
  # print(noquote(c("Twod: ", twod)));
  
  return (twod) 
}



second_term <- function(alpha, beta, m0, m1)
{
#
# this function computes the second term of the three terms
# in the expression for the variance:
#
  oned_2_integrant <- function(x1, x2, alpha, beta, m0, m1) 
  {
    
    temp = (x1 * x2 * alpha  + 1 - alpha)^(m0 - 2) *
      (x1 * x2 * beta + 1 - beta)^m1;
    
    temp = temp * x1 * x2;
    
    return (temp);
  }
  
  
  
  oned_2_integral <- function(x2, alpha, beta, m0, m1)
  {
    res = integrate(oned_2_integrant, lower=0, upper= 0.6, x2 = x2, alpha=alpha, beta = beta, m0 = m0, m1 = m1)$value
    
    res = res + integrate(oned_2_integrant, lower=0.6, upper= 0.7, x2 = x2, alpha=alpha, beta = beta, m0 = m0, m1 = m1)$value
    
    res = res + integrate(oned_2_integrant, lower=0.7, upper= 0.8, x2 = x2, alpha=alpha, beta = beta, m0 = m0, m1 = m1)$value
    
    res = res + integrate(oned_2_integrant, lower=0.8, upper= 0.9, x2 = x2, alpha=alpha, beta = beta, m0 = m0, m1 = m1)$value
    
    res = res + integrate(oned_2_integrant, lower=0.9, upper= 0.91, x2 = x2, alpha=alpha, beta = beta, m0 = m0, m1 = m1)$value
    
    res = res + integrate(oned_2_integrant, lower=0.91, upper= 0.95, x2 = x2, alpha=alpha, beta = beta, m0 = m0, m1 = m1)$value
    
    res = res + integrate(oned_2_integrant, lower=0.95, upper= 0.96, x2 = x2, alpha=alpha, beta = beta, m0 = m0, m1 = m1)$value
    
    res = res + integrate(oned_2_integrant, lower=0.96, upper= 0.97, x2 = x2, alpha=alpha, beta = beta, m0 = m0, m1 = m1)$value
    
    res = res + integrate(oned_2_integrant, lower=0.97, upper= 0.98, x2 = x2, alpha=alpha, beta = beta, m0 = m0, m1 = m1)$value
    
    res = res + integrate(oned_2_integrant, lower=0.98, upper= 0.99, x2 = x2, alpha=alpha, beta = beta, m0 = m0, m1 = m1)$value
    
    res = res + integrate(oned_2_integrant, lower=0.99, upper= 0.995, x2 = x2, alpha=alpha, beta = beta, m0 = m0, m1 = m1)$value
    
    #res = res + integrate(oned_2_integrant, lower=0.995, upper= 1.0, x2 = x2, alpha=alpha, beta = beta, m0 = m0, m1 = m1)$value
    res = res + integrate(oned_2_integrant, lower=0.995, upper= 0.9998, x2 = x2, alpha=alpha, beta = beta, m0 = m0, m1 = m1, subdivisions = 10000)$value
    res = res + integrate(oned_2_integrant, lower=0.9998, upper= 1.0, x2 = x2, alpha=alpha, beta = beta, m0 = m0, m1 = m1, subdivisions = 10000)$value
    
    return (res)
  }
  
  twod_2 = integrate(Vectorize(oned_2_integral), alpha=alpha, beta=beta, m0=m0, m1=m1, lower=0, upper=0.3)$value
  
  twod_2 = twod_2 + integrate(Vectorize(oned_2_integral), alpha=alpha, beta=beta, m0=m0, m1=m1, lower=0.3, upper=0.4)$value
  
  twod_2 = twod_2 + integrate(Vectorize(oned_2_integral), alpha=alpha, beta=beta, m0=m0, m1=m1, lower=0.4, upper=0.5)$value
  
  twod_2 = twod_2 + integrate(Vectorize(oned_2_integral), alpha=alpha, beta=beta, m0=m0, m1=m1, lower=0.5, upper=0.6)$value
  
  twod_2 = twod_2 + integrate(Vectorize(oned_2_integral), alpha=alpha, beta=beta, m0=m0, m1=m1, lower=0.6, upper=0.7)$value
  
  twod_2 = twod_2 + integrate(Vectorize(oned_2_integral), alpha=alpha, beta=beta, m0=m0, m1=m1, lower=0.7, upper=0.8)$value
  
  twod_2 = twod_2 + integrate(Vectorize(oned_2_integral), alpha=alpha, beta=beta, m0=m0, m1=m1, lower=0.8, upper=0.9)$value
  
  twod_2 = twod_2 + integrate(Vectorize(oned_2_integral), alpha=alpha, beta=beta, m0=m0, m1=m1, lower=0.9,  upper=0.925)$value
  
  twod_2 = twod_2 + integrate(Vectorize(oned_2_integral), alpha=alpha, beta=beta, m0=m0, m1=m1, lower=0.925, upper=0.95)$value
  
  twod_2 = twod_2 + integrate(Vectorize(oned_2_integral), alpha=alpha, beta=beta, m0=m0, m1=m1, lower=0.95, upper=0.96)$value
  
  twod_2 = twod_2 + integrate(Vectorize(oned_2_integral), alpha=alpha, beta=beta, m0=m0, m1=m1, lower=0.96,  upper=0.97)$value
  
  twod_2 = twod_2 + integrate(Vectorize(oned_2_integral), alpha=alpha, beta=beta, m0=m0, m1=m1, lower=0.97,  upper=0.98)$value
  
  twod_2 = twod_2 + integrate(Vectorize(oned_2_integral), alpha=alpha, beta=beta, m0=m0, m1=m1, lower=0.98,  upper=0.99)$value
  
  #twod_2 = twod_2 + integrate(Vectorize(oned_2_integral), alpha=alpha, beta=beta, m0=m0, m1=m1, lower=0.99, upper=1.)$value
  
  twod_2 = twod_2 + integrate(Vectorize(oned_2_integral), alpha=alpha, beta=beta, m0=m0, m1=m1, lower=0.99,  upper=0.9998)$value
  
  twod_2 = twod_2 + integrate(Vectorize(oned_2_integral), alpha=alpha, beta=beta, m0=m0, m1=m1, lower=0.9998, upper=1.)$value
  
  
  # print(c(integrate(Vectorize(oned_2_integral), alpha=alpha, beta=beta, m0=m0, m1=m1, lower=0, upper=1, subdivisions = 100000)$value,
  # twod_2));
  
  twod_2 = twod_2 * m0 * (m0- 1) * alpha^2  / (1 - (1 - alpha)^m0*(1 - beta)^m1);
  
  #print(noquote(c("Twod_2: ", twod_2)));
  
  return (twod_2) 
  
}


third_term <- function(alpha, beta, m0, m1)
{

#
# this function computes the third term of the three terms
# in the expression for the variance. It is the square of
# the pFDR. It is a one-demensional integral
#
# The script first determines the first value of the argument (x) for
# which the integrand is larger than 10^(-7). Then it does the integration
# for two ranges: up to the determined value and after that.
#
#

  integrant <- function(s, alpha, beta, m0, m1) 
  {
    (s * alpha  + 1 - alpha)^(m0 - 1) *
      (s * beta + 1 - beta)^m1 ;
  }
  
  
  pFDR_theoretical <- integrate(integrant, lower=0, upper = 0.1, alpha = alpha, beta = beta, m0 = m0, m1 = m1)$value +
    integrate(integrant, lower=0.1, upper = .20, alpha = alpha, beta = beta, m0 = m0, m1 = m1)$value +
    integrate(integrant, lower=0.2, upper = .30, alpha = alpha, beta = beta, m0 = m0, m1 = m1)$value +
    integrate(integrant, lower=0.3, upper = .40, alpha = alpha, beta = beta, m0 = m0, m1 = m1)$value + 
    integrate(integrant, lower=0.4, upper = .50, alpha = alpha, beta = beta, m0 = m0, m1 = m1)$value +
    integrate(integrant, lower=0.5, upper = .60, alpha = alpha, beta = beta, m0 = m0, m1 = m1)$value +
    integrate(integrant, lower=0.6, upper = .70, alpha = alpha, beta = beta, m0 = m0, m1 = m1)$value +
    integrate(integrant, lower=0.70, upper = .80, alpha = alpha, beta = beta, m0 = m0, m1 = m1)$value +
    integrate(integrant, lower=0.80, upper = .90, alpha = alpha, beta = beta, m0 = m0, m1 = m1)$value +
    integrate(integrant, lower=0.9, upper = .95, alpha = alpha, beta = beta, m0 = m0, m1 = m1)$value +
    integrate(integrant, lower=0.95, upper = .96, alpha = alpha, beta = beta, m0 = m0, m1 = m1)$value +
    integrate(integrant, lower=0.96, upper = .97, alpha = alpha, beta = beta, m0 = m0, m1 = m1)$value +
    integrate(integrant, lower=0.97, upper = .98, alpha = alpha, beta = beta, m0 = m0, m1 = m1)$value +
    integrate(integrant, lower=0.98, upper = .99, alpha = alpha, beta = beta, m0 = m0, m1 = m1)$value +
    #integrate(integrant, lower=0.99, upper = 1., alpha = alpha, beta = beta, m0 = m0, m1 = m1)$value 
    integrate(integrant, lower=0.99, upper = 0.9998, alpha = alpha, beta = beta, m0 = m0, m1 = m1, subdivisions=10000)$value +
    integrate(integrant, lower=0.9998, upper = 1., alpha = alpha, beta = beta, m0 = m0, m1 = m1)$value
  
  
  split_value = return_limit(alpha, beta, m0, m1);
  
  #print(c("pFDR before: ", pFDR_theoretical));
  
  if(split_value > 0.0)
  {
    pFDR_theoretical <- integrate(integrant, lower=0, upper = split_value, alpha = alpha, beta = beta, m0 = m0, m1 = m1)$value +
      integrate(integrant, lower=split_value, upper = 1., alpha = alpha, beta = beta, m0 = m0, m1 = m1)$value;
    
  }
  
  #print(c("pFDR after: ", pFDR_theoretical));
  
  
  #   print("Integral only: ");
  #   print(pFDR_theoretical );
  #   print("Integral X m0 x alpha");
  #   print(pFDR_theoretical * m0 * alpha );
  
  
  pFDR_theoretical <-  pFDR_theoretical * m0 * alpha / (1. - (1. - alpha)^m0 * (1. - beta)^m1 );
  
  #
  #  for the variance, the square is needed
  #
  #   print(pFDR_theoretical^2);
  return (pFDR_theoretical^2)
}



pFDR_variance <- function(alpha, beta, m0, m1)
{
#
# the variance is the made of the three terms
#
#
  
  temp3 = third_term(alpha, beta, m0, m1);
  
  temp = first_term(alpha, beta, m0, m1);
  
  temp = temp + second_term(alpha, beta, m0, m1);
  
  temp = temp - temp3;
  
  if(temp < 0.0)
  {
    print("The integral calculations went wrong: Need to increase the grid density");
    
    return (-1);
  }
  
  options(digits= 5);
  #format(print(noquote(c("pFDR, variance,SD: ", sqrt(temp3), temp, sqrt(temp) ) ), digits=4));
  
  return(temp);
  
}


#
#   Calculations using integral2 function
#   from pracma package.
#
#

# 
# library(pracma)
# 
# for_oned_integral <- function(x, alpha, beta, m0, m1)
# {
#   
#   # alpha = 0.005; beta = 0.6764718; m0 = 60000; m1 = 40000;
#   
#   
#   (x * alpha  + 1 - alpha)^(m0 - 1) *
#     (x * beta + 1 - beta)^m1 ;
#   
# }
# 
# 
# for_integral2 <- function(x, y, alpha, beta, m0, m1)
# {
#   
#   # alpha = 0.005; beta = 0.6764718; m0 = 60000; m1 = 40000;
#   
#   
#   x * y * (x * y * alpha  + 1 - alpha)^(m0 - 2) *
#     (x * y * beta + 1 - beta)^m1;
#   
# }
# 
# 
# for_integral2_second <- function(x, y, alpha, beta, m0, m1)
# {
#   
#   #alpha = 0.005; beta = 0.6764718; m0 = 60000; m1 = 40000;
#   
#   
#   (x * y * alpha + 1 - alpha)^(m0-1) * (x* y * beta + 1 - beta)^m1;
#   
# }
# 
# temp1 = integral2(for_integral2, 0, 1, 0, 1, alpha = alpha, beta = beta, m0 = m0, m1 = m1, abstol = 0.0, reltol = 1e-11, maxlist = 50000)$Q * m0 * (m0-1)*alpha^2
# 
# temp1 = temp1 / (1 - (1 - alpha)^m0 * (1 - beta)^m1);
# 
# temp2 = integral2(for_integral2_second, 0, 1, 0, 1, alpha = alpha, beta = beta, m0 = m0, m1 = m1, abstol = 0.0, reltol = 1e-11, maxlist = 50000)$Q * m0 *alpha
# 
# temp2 = temp2 / (1 - (1 - alpha)^m0 * (1 - beta)^m1);
# 
# temp3 = integral(for_oned_integral, 0, 1, alpha = alpha, beta = beta, m0 = m0, m1 = m1, abstol = 0.0, reltol = 1e-11, no_intervals = 1000) * m0 *alpha
# 
# temp3 = temp3 / (1 - (1 - alpha)^m0 * (1 - beta)^m1);
# 
# variance = temp1 + temp2 - temp3^2