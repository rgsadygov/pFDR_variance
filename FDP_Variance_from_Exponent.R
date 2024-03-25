#
#
# An R script to compute the false discovery rate and the variance of the
# of the false discovery proportions using the results from peptide
# identifications with forward and reversed sequences using the approximate
# formulas in "Exact Integral Formulas for False Discovery Rate
#   and the Variance of False Discovery Proportion"
#
# After sourcing this is R code, the function is to be called as:
#
# pFDR_Exponent_variance(N_falsePositives, N_significant, N_sample) 
#
# For example, if
# N_falsePositves = 27294; N_significant = 107799; N_sample = 909800 + 178900 
#
# Then:
# pFDR_Exponent_variance(27294, 107799, 1088700)
#
# In the output: 2.5319e-01 1.3192e-03 1.7402e-06
# FDR = 2.5319e-01; standard deviation of FDP = 1.3192e-03
# variance of FDP = 1.7402e-06
#
#

N_grid    = 10;

oned_integrant_exponential <- function(x1, x2, N_significant)
{
	 temp = exp(-N_significant*(1 - x1 * x2));

	 return (temp);
}


oned_2_integrant <- function(x1, x2, N_significant) 
{
         temp = exp(- N_significant * (1 - x1*x2));

	 temp = temp * x1 * x2;

	 return (temp);
}

#
#   This function determines the value of the integration variable at which the exponential
#   function becomes non-zero.
#   Its arguments are: N_significant, and the number of sample points (N_crit);
#   N_significant - the number of trials called significant (passing a certain threshold)
#
critical_point_1d <- function(N_significant, N_crit)
{  
    for(i in 1:N_crit)
    {
       delta_x = i / N_crit;

       temp = exp(-N_significant*(1 - delta_x) );

       if(temp > 10^(-7) )
       {
              break;
       }
    }

   mid_step = i / N_crit;

#   print(c(i, mid_step, N_significant));   

   return(mid_step)
   
}  


#
#   The function computes the first part of the first integral in Eq. (13) of the
#   supporting note in "Exact Integral Formulas for False Discovery Rate
#   and the Variance of False Discovery Proportion"
#
#   Its arguments are: N_false, N_significant, and mid_step 
#   N_false - the number of false positives
#   N_significant - the number of trials called significant (passing a certain threshold)
#   mid_step the integration point at which the integrant is non-zero.
#
first_term_exponential <- function(N_false, N_significant, mid_step)
{

    delta = (1. - mid_step)/N_grid;

    low_point = mid_step; 

    oned_integral <- function(x2, N_significant)
    {

	res = integrate(oned_integrant_exponential, lower=0, upper= mid_step, x2 = x2, N_significant = N_significant)$value

        delta = (1. - mid_step)/N_grid;

        low_point = mid_step; 

        for(i in 1:N_grid)
	{

	   upper_point = mid_step + i * delta;

          # print(c(low_point, upper_point));
	   
	   res = res + integrate(oned_integrant_exponential, lower=low_point, upper= upper_point, x2 = x2, N_significant = N_significant)$value

	   low_point = upper_point;
	}

        return (res)
    }
 

    delta = (1. - mid_step)/N_grid;

    low_point = mid_step; 
    
    twod = integrate(Vectorize(oned_integral), N_significant=N_significant, lower=0, upper=mid_step)$value

    for(i in 1:N_grid)
    {
        upper_point = mid_step + i * delta;
	
	twod = twod + integrate(Vectorize(oned_integral), N_significant=N_significant, lower=low_point, upper=upper_point)$value

	low_point = upper_point;
    }

    twod = twod * N_false;

    
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

second_term_exponential <- function(N_false, N_significant, mid_step)
{

    delta = (1. - mid_step) / N_grid;

    low_point = mid_step; 

    oned_2_integral <- function(x2, N_significant)
    {
        res = integrate(oned_2_integrant, lower=0, upper= mid_step, x2 = x2, N_significant = N_significant)$value

	low_point = mid_step; 

        for(i in 1:N_grid)
        {
           upper_point = mid_step + i * delta;
	
           res = res + integrate(oned_2_integrant, lower=low_point, upper= upper_point, x2 = x2, N_significant = N_significant)$value

	   low_point = upper_point;
        }


        return (res)
    }

    twod_2 = integrate(Vectorize(oned_2_integral), N_significant = N_significant, lower=0, upper=mid_step)$value

    low_point = mid_step; 

    for(i in 1:N_grid)
    {
       upper_point = mid_step + i * delta;
	
       twod_2 = twod_2 + integrate(Vectorize(oned_2_integral), N_significant = N_significant, lower=low_point, upper=upper_point)$value

       low_point = upper_point;
    }

    twod_2 = twod_2 * N_false^2;

    return (twod_2) 

}


#
# The function computes the FDP variance in Eq. (13) of the
#   supporting note in "Exact Integral Formulas for False Discovery Rate
#   and the Variance of False Discovery Proportion" 
#   N_sample - the number of all (null and alternative) hypothesis tests
#   N_false - the number of null hypotheses call significant
#   N_signifcant - the number of all (null and alternative) hypotheses called significant
#

pFDR_Exponent_variance <- function(N_false, N_significant, N_sample) 
{

    mid_step = 1;

    N_crit = N_sample;

    while (mid_step == 1)
    {
        mid_step = critical_point_1d(N_significant, N_crit);

	# print(mid_step);

	N_crit = 10 * N_crit;
    }

   if(mid_step > 10)
   {
       mid_step = mid_step - 10;
   }

   temp = first_term_exponential(N_false, N_significant, mid_step);

   # print(c("First part from Exponent ", temp));

   temp = temp + second_term_exponential(N_false, N_significant, mid_step);

   # print(c("Second part from Exponent ", second_term_exponential(N_false, N_significant, mid_step)) );

   # print(N_false / N_significant);

   temp = temp - (N_false / N_significant) ^2;

   # print(c("Third part from Exponent ", (N_false / N_significant) ^2) );

   # print("The FDP, SD(FDP), and var(FDP), from exponential approximation, are: ");

   # print(c((N_false/N_significant), sqrt(temp), temp));
   
   return(temp);
 
}

