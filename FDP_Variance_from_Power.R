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
#

third_term_power <- function(alpha, beta, m0, m1, mid_step)
{
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
   return (pFDR_theoretical^2)
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

   temp = first_term_power(alpha, beta, m0, m1, mid_step);

   # print(c("First part from power ", temp));

   temp = temp + second_term_power(alpha, beta, m0, m1, mid_step);

   # print(c("Second part from power ", second_term_power(alpha, beta, m0, m1, mid_step)) );

   third_term = third_term_power(alpha, beta, m0, m1, mid_step);

   temp = temp - third_term;

   # print(c("Third part from power ", third_term ));

   # print("The FDP, SD(FDP), and var(FDP), from exponential approximation, are: ");

   # print(c(sqrt(third_term), sqrt(temp), temp));
   
   return (temp);
   
}

