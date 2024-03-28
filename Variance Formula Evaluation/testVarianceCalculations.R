library("ggplot2")


##=======================================================
##======code for computing variance of FDR===============
##=======================================================

#
#
# An R script to compute the false discovery rate and the variance
# of the false discovery proportions using
# probabilistic properties, such as the type I (alpha) and
# type II (1- beta) errors, and the numbers of null (m0)
# and alternative hypotheses (m1).
# After sourcing this R code, the function is to be called as:
# pFDR_Power_variance(alpha, beta, mo, m1) 
#
# For example if
# alpha = 0.03, beta= 0.45, m0 = 909800, m1 = 178900
#
# Then:
# pFDR_Power_variance(0.03, 0.45, 909800, 178900)
# In the output: 2.5319e-01 1.2591e-03 1.5854e-06
# FDR = 2.5319e-01; standard deviation of FDP = 1.2591e-03
# variance of FDP = 1.5854e-06
#

N_grid    = 10;

oned_integrant_power <- function(x1, x2, alpha, beta, m0, m1)
{
	 temp = (x1 * x2 * alpha  + 1 - alpha)^(m0 - 1) *
         (x1 * x2 * beta + 1 - beta)^m1

	 return (temp);
}


oned_2_integrant_power <- function(x1, x2, alpha, beta, m0, m1) 
{

	 temp = (x1 * x2 * alpha  + 1 - alpha)^(m0 - 2) *
         (x1 * x2 * beta + 1 - beta)^m1;

	 temp = temp * x1 * x2;

	 return (temp);
}

#
#   This function determines the value of the integration variable at which the exponential
#   function becomes non-zero.
#   alpha, (1- beta) - are the probabilites of type I and II errros.
#   m0 and m1 are the numbers of the null and alternative hypotheses
#   N_crit - is the number of sample points. It is increased to find an
#   appropriate mid_step - is the value of the variable in the integrand
#   at which the integrand is non-zero ( > 1 ^(-7)).
#
#
critical_point_1d <- function(alpha, beta, m0, m1, N_crit)
{  
    for(i in 1:N_crit)
    {
       delta_x = i / N_crit;

       temp = (delta_x* alpha  + 1 - alpha)^(m0 - 1) *
         (delta_x * beta + 1 - beta)^m1

       if(temp > 10^(-7) )
       {
              break;
       }
    }

   mid_step = i / N_crit;

   return(mid_step)
   
}  


#
#   The function computes the first part of the first integral in Eq. (12) of the
#   supporting note in "Exact Integral Formulas for False Discovery Rate
#   and the Variance of False Discovery Proportion"
#
#   alpha, (1- beta) - are the type I and II errros.
#   m0 and m1 are the numbers of the null and alternative hypotheses
#   mid_step - is the value of the variable in the integrand
#   at which the integrand is non-zero.
#
first_term_power <- function(alpha, beta, m0, m1, mid_step)
{

    delta = (1. - mid_step)/N_grid;

    low_point = mid_step; 

    oned_integral <- function(x2, alpha, beta, m0, m1)
    {

	res = integrate(oned_integrant_power, lower=0, upper= mid_step, x2 = x2, alpha=alpha, beta = beta, m0 = m0, m1 = m1)$value

        delta = (1. - mid_step)/N_grid;

        low_point = mid_step; 

        for(i in 1:N_grid)
	{

	   upper_point = mid_step + i * delta;

	   res = res + integrate(oned_integrant_power, lower=low_point, upper= upper_point, x2 = x2, alpha=alpha, beta = beta, m0 = m0, m1 = m1)$value; 

	   low_point = upper_point;
	}

        return (res)
    }
 

    delta = (1. - mid_step)/N_grid;

    low_point = mid_step; 
    
    twod = integrate(Vectorize(oned_integral), alpha=alpha, beta=beta, m0=m0, m1 = m1, lower=0, upper=mid_step)$value

    for(i in 1:N_grid)
    {
        upper_point = mid_step + i * delta;
	
        twod = twod + integrate(Vectorize(oned_integral), alpha=alpha, beta=beta, m0=m0, m1=m1, lower=low_point, upper=upper_point)$value

	low_point = upper_point;
    }

    twod = twod * m0 * alpha  / (1 - (1 - alpha)^m0*(1 - beta)^m1);

    return (twod) 
    
}

#
#
#   The function computes the second integral in Eq. (S13) of the
#   supporting note in "Exact Integral Formulas for False Discovery Rate
#   and the Variance of False Discovery Proportion" 
#   alpha, (1- beta) - are the type I and II errros.
#   m0 and m1 are the numbers of the null and alternative hypotheses
#   mid_step - is the value of the variable in the integrand
#   at which the integrand is non-zero.
#

second_term_power <- function(alpha, beta, m0, m1, mid_step)
{

    delta = (1. - mid_step) / N_grid;

    low_point = mid_step; 

    oned_2_integral_power <- function(x2, alpha, beta, m0, m1)
    {
        res = integrate(oned_2_integrant_power, lower=0, upper= mid_step, x2 = x2, alpha=alpha, beta=beta, m0=m0, m1=m1)$value

	low_point = mid_step; 

        for(i in 1:N_grid)
        {
           upper_point = mid_step + i * delta;
	
           res = res + integrate(oned_2_integrant_power, lower=low_point, upper= upper_point, x2 = x2, alpha=alpha, beta=beta, m0=m0, m1=m1)$value

	   low_point = upper_point;
        }


        return (res)
    }

    twod_2 = integrate(Vectorize(oned_2_integral_power), alpha=alpha, beta=beta, m0=m0, m1=m1, lower=0, upper=mid_step)$value

    low_point = mid_step; 

    for(i in 1:N_grid)
    {
       upper_point = mid_step + i * delta;
	
       twod_2 = twod_2 + integrate(Vectorize(oned_2_integral_power), alpha=alpha, beta=beta, m0=m0, m1=m1, lower=low_point, upper=upper_point)$value

       low_point = upper_point;
    }

    twod_2 = twod_2 * m0 * (m0-1) * alpha^2 / (1 - (1 - alpha)^m0*(1 - beta)^m1);

    return (twod_2) 

}




#
# this function computes the third term of the three terms
# in the expression for the variance. It is the square of
# the pFDR. It is a one-demensional integral
#third_term_power <- function(alpha, beta, m0, m1, mid_step)

pFDR_IntegralApproach <- function(alpha, beta, m0, m1, mid_step= -Inf)
{
   if (mid_step == -Inf){
     mid_step = get_mid_step(alpha, beta, m0, m1);
   }
   integrant <- function(s, alpha, beta, m0, m1) 
    {
         (s * alpha  + 1 - alpha)^(m0 - 1) *
         (s * beta + 1 - beta)^m1 ;
    }

    delta = (1. - mid_step) / N_grid;

    pFDR_theoretical <- integrate(integrant, lower=0, upper = mid_step, alpha = alpha, beta = beta, m0 = m0, m1 = m1)$value;

    low_point = mid_step; 

    for(i in 1:N_grid)
    {
       upper_point = mid_step + i * delta;
	
       pFDR_theoretical = pFDR_theoretical + integrate(integrant, lower=low_point, upper = upper_point, alpha = alpha, beta = beta, m0 = m0, m1 = m1)$value

       low_point = upper_point;
    }


   pFDR_theoretical <-  pFDR_theoretical * m0 * alpha / (1. - (1. - alpha)^m0 * (1. - beta)^m1 );

#
#  for the variance, the square is needed
#
   return (pFDR_theoretical)
}


get_mid_step <- function(alpha, beta, m0, m1){
	mid_step = 1;

    N_crit = m0 + m1;

    while (mid_step == 1)
    {
        mid_step = critical_point_1d(alpha, beta, m0, m1, N_crit);

	# print(mid_step);

	N_crit = 10 * N_crit;
    }

   if(mid_step > 10)
   {
       mid_step = mid_step - 10;
   }
   return (mid_step);
}
#
# The function computes the FDP variance in Eq. (12) of the
#   supporting note in "Exact Integral Formulas for False Discovery Rate
#   and the Variance of False Discovery Proportion"
#   alpha, (1- beta) - are the type I and II errros.
#   m0 and m1 are the numbers of the null and alternative hypotheses
#
#

pFDR_Power_variance <- function(alpha, beta, m0, m1)
{

   mid_step = get_mid_step(alpha, beta, m0, m1);

   temp = first_term_power(alpha, beta, m0, m1, mid_step);

   # print(c("First part from power ", temp));

   temp = temp + second_term_power(alpha, beta, m0, m1, mid_step);

   # print(c("Second part from power ", second_term_power(alpha, beta, m0, m1, mid_step)) );

   third_term = pFDR_IntegralApproach(alpha, beta, m0, m1, mid_step);

   temp = temp - third_term^2;

   # print(c("Third part from power ", third_term ));

   # print("The FDP, SD(FDP), and var(FDP), from exponential approximation, are: ");

   # print(c((third_term), sqrt(temp), temp));
   
   return (temp);
   
}

##=======================================================
##=============End var pFDR==============================
##=======================================================



##=======================================================
##=========Evaluation of variance calculation============
##=======================================================

#determine the combined density at a given n choose k (null distribution) and n choose j (alternative distribution)
singleRoundDensity <- function(alpha, beta, m0, m1, k, j) {
  nullDensity = dbinom(k, m0, alpha) #density for null Distribution
  altDensity = dbinom(j, m1, beta) #density for alternative Distribution
  combinedDensity = nullDensity * altDensity #combined Density
  return(combinedDensity)
}

#calculate the variance of the FDP given a set of parameters 
calculateVarQSummation <- function(alpha, beta, m0, m1) {
  alphaBase = (1-alpha)^m0
  betaBase = (1-beta)^m1
  normalization = 1/(1-alphaBase*betaBase) #calculates the normalization value used in the Variance of FDP formula
  
  expectationOfSqr = 0 #running sum of the expectation of the square of the FDP
  sqrOfExpectation = 0 #running sum of the square of the expectation of the FDP
  
  for(k in 1:m0) { #outer loop for summation over k (null distribution sample size)
    for (j in 1:m1) { #inner loop for summation over j (alternative distribution sample size)
      indexFaction = k/(k+j) #calculate ratio of k to k+j
      
      density = singleRoundDensity(alpha, beta, m0, m1, k, j) #retrieve combined density
      expectationOfSqr = expectationOfSqr + (indexFaction^2 * density) #update running sum of expectation of square
      sqrOfExpectation = sqrOfExpectation + (indexFaction * density) #update running sum of square of expectation
    }
  }
  
  qValueVar = (normalization * expectationOfSqr) - (normalization * sqrOfExpectation)^2 #calculate Final Value
  return(qValueVar)
} 

#loops through ranges on different parameters and calculatse the relative difference between 
#the FDP variance calculated from binomial summation and integration approach, respective to the the binomial value.
iterateOverParams <- function() {
  m0_range = seq(50, 100, 10) #different paramter ranges
  m1_range = seq(25, 50, 5) 
  alpha_range = seq(0.01, 0.05, 0.01) 
  beta_range = seq(0.4, 0.8, 0.1) 
  
  columnNames <- c("alpha", "beta", "m0", "m1",  "VarFDP_Binomal", "VarFDP_Integration")
  toleranceFrame <- data.frame(matrix(ncol=length(columnNames), nrow=0))
  colnames(toleranceFrame) <- columnNames #generate dataframe to store parameters and calculate values
  
  for(alpha in alpha_range) {
    for(beta in beta_range) {
      for(m0 in m0_range) {
        for(m1 in m1_range) {
          
          VarQSummation <- calculateVarQSummation(alpha, beta, m0, m1) #Variance of FDP from binomial approach
          VarQIntegral <- pFDR_Power_variance(alpha, beta, m0, m1) #Variance of FDP from integration approach
          
          newRow = c(alpha, beta, m0, m1, VarQSummation, VarQIntegral)
          toleranceFrame[nrow(toleranceFrame) + 1,] = newRow #update dataframe
        }
      }
    }
  }
  
  print("finished")
  return(toleranceFrame)
}


#dumps results to csv and generates figure

#dumps results to csv and generates figure
plotLinearGraphAndDumpOutput <- function(toleranceFrame) {
  #outFile <- "C:/Users/Justin/Desktop/FDR_pFDR/VarFDP_BinomalandIntegration_Differences_New.csv"
  #write.csv(toleranceFrame, outFile, row.names=FALSE) #dump results to csv file
  
  #generate graph plotting Variance FDP from binomial approach (x) against integration approach (y)
  #captionLine <- "Comparision of the FDP Variance Values Calculated from the Binomial Approach (x-axis) and Integration Approach (y-axis). Red Line Signifies Equality (x=y)"
  
  ggplot(toleranceFrame, aes(x=VarFDP_Binomal, y=VarFDP_Integration)) + 
    geom_point(size=2) +
    geom_abline(slope=1, intercept=0, 
                col="red", linewidth=1) +
    
    xlab("\nVariance of FDP Calculated \nFrom Binomial Formula") +
    ylab("Variance of FDP Calculated \nFrom Integration Formula\n") +
    
    theme_bw() + #remove the grey background and grid lines
    theme(panel.border = element_rect(linewidth = 1),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.title = element_text(color = "black",
                                    size = 16),
          axis.text = element_text(color = "black",
                                   size = 16)) 
}


toleranceFrame <- iterateOverParams()
plotLinearGraphAndDumpOutput(toleranceFrame)

##=======================================================
##================End Comparision========================
##=======================================================