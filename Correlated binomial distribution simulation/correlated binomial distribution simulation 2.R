#======================================================================
#======================================================================
#  FDP computation from  correlated binomial random variables
#  Source: Kotz, Samuel, et al. Encyclopedia of Statistical Sciences, Volume 1. John Wiley & Sons, 2005.
#         page no. 582

#  The code below implements FDP computation from correlated binominal distribution.
#  Suppose X (corresponds to the null)  is a Bernoulli(P(X=1) = alpha) variable (P(X = 0) = 1- alpha), and
#  Y (corresponds to the alternative) is another Bernoulli(P(Y=1)= beta) variable (P(Y=0) = 1 - beta).
#  To generate their joint distribution we specified
#  all four combinations of outcomes as follows (similar to the above mentioned book).
#     P(X=0,Y=0) = a,                       # denote this combination as "(0,0)"
#     P(X=1,Y=0) = 1- beta -a,              # denote this combination as "(1,0)"
#     P(X=0,Y=1) = 1- alpha - a,            # denote this combination as "(0,1)"
#     P(X=1,Y=1) = a + alpha + beta -1,     # denote this combination as "(1,1)"
#
# Where "a" is computed from the correlation coefficient (rho) as follows:
#   a = (1- alpha )(1- beta) + rho *  sqrt (alpha * beta *(1- alpha )(1- beta)  )
#
#
# The FDP is computed from simulation (using sample() in R) by  sampling n=10000 points using the
# probabilies of the four outcomes: (P(X = 1, Y = 1), P(X = 1, Y = 0),P(X = 0, Y = 1), P(X = 0, Y = 0)).
#
#   Then:
#
#   FP =  number of outcomes of "(1,0)" + number of outcomes of "(1,1)"
#   TP =  number of outcomes of "(0,1)" + number of outcomes of "(1,1)"
#   FDP = FP / (FP + TP)

#======================================================================
#======================================================================
#======================================================================



#======================================================================
#========================FDP variance==================================
#======================================================================
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

#======================================================================
#====================END FDP variance==================================
#======================================================================


#
# The function, compute_fdp) computes the false discovery proportion given, rho(correlation),
# alpha (significane level of the null), and beta (the power of the test)
#
#
compute_fdp = function (rho, alpha, beta, n.sim) {
  # Compute p(x=0, y=0) from the correlation value
  a.0 <-
    rho * sqrt(alpha * beta * (1 - alpha) * (1 - beta)) + (1 - alpha) * (1 -
                                                                           beta)
  
  # Compute the four probabilities for the joint distribution.
  prob <-
    c(
      `(0,0)` = a.0,
      `(1,0)` = 1 - beta - a.0,
      `(0,1)` = 1 - alpha - a.0,
      `(1,1)` = a.0 + alpha + beta - 1
    )
  if (min(prob) < 0) {
    print(prob)
    stop("Error: a probability is negative.")
  }
  
  #
  # generation of correlated Binomial variables.
  #
  
  outcomes <-
    sample.int(4, n.sim , replace = TRUE, prob = prob) #This will generate random numbers between 1 and 4 with the
  #probabilities assigned above. where
  #     1 denote the combination  "(0,0)"
  #     2 denote the combination  "(1,0)"
  #     3 denote the combination  "(0,1)"
  #     4 denote the combination  "(1,1)"
  
  
  
  # get the numbers of P[X = 1, Y = 1], P[X = 1, Y = 0],P[X = 0, Y = 1], and = P[X = 0, Y = 0]
  p_x0_y0 = 0   # p(x=0,y=0)
  p_x0_y1 = 0   # p(x=0,y=1)
  p_x1_y0 = 0   # p(x=1,y=0)
  p_x1_y1 = 0   # p(x=1,y=1)
  
  for (i  in seq(1:length(outcomes))) {
    if (outcomes[i] == 1) {
      p_x0_y0 = p_x0_y0 + 1
    }
    else if (outcomes[i] == 3) {
      p_x0_y1 = p_x0_y1 + 1
    }
    else if (outcomes[i] == 2) {
      p_x1_y0 = p_x1_y0 + 1
    }
    else if (outcomes[i] == 4) {
      p_x1_y1 = p_x1_y1 + 1
    }
    
  }
  
  # compute fdp
  FP =  (p_x1_y1 + p_x1_y0)    # number of false positives (FP).
  TP =  (p_x1_y1 + p_x0_y1)     # number of true positives (TP).
  fdp = FP / (FP + TP)
  
  # print(c('rho ', rho, 'FDP ', fdp))
  return(fdp)
}


#======================================================================

# run the simulation for diffrent Correlation values
rho = c(-0.99, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.2, 0.1, -0.0001)
alpha = 0.1
beta = 0.9
n.sim <- 10000

all_res = list()
for (iter in seq(1, 100)) {
  print(iter)
  fdps = c(rep(-1, length(rho)))
  for (i in seq(1, length(rho))) {
    fdps[i] = compute_fdp(rho[i], alpha, beta, n.sim)
  }
  all_res = append(all_res, fdps)
}


# compute the integeral FDR
c_fdr = pFDR_IntegralApproach (alpha, beta, n.sim, n.sim)
c_var = pFDR_Power_variance(alpha, beta, n.sim, n.sim)

# prepare data for plot
df = data.frame(
  unlist(all_res[seq(1, 1000, 10)]),
  unlist(all_res[seq(2, 1000, 10)]),
  unlist(all_res[seq(3, 1000, 10)]),
  unlist(all_res[seq(4, 1000, 10)]),
  unlist(all_res[seq(5, 1000, 10)]),
  unlist(all_res[seq(6, 1000, 10)]),
  unlist(all_res[seq(7, 1000, 10)]),
  unlist(all_res[seq(8, 1000, 10)]),
  unlist(all_res[seq(9, 1000, 10)]),
  unlist(all_res[seq(10, 1000, 10)])
)

colnames(df) = c(rho[1], rho[2], rho[3],
                 rho[4], rho[5], rho[6],
                 rho[7], rho[8], rho[9], rho[10])


df2 = data.frame(rho, rep(c_fdr, length(rho)), rep(c_var, length(rho)))
colnames(df2) = c('rho', 'FDR', 'VAR_FDR')


#======================================================================
# Plot result


plot_errorbars <- function(x, y, err, ylim = NULL) {
  arrows(
    x,
    y - err,
    x,
    y + err,
    length = 0.05,
    angle = 90,
    code = 3,
    col = 'red'
  )
}


boxplot(
  df,
  col = 'white',
  ylab = 'FDP',
  xlab = 'Correlation coefficient',
  ylim = c(min(fdps) - 0.005, max(fdps) + 0.005),
  cex = 1.5,
  cex.axis = 1.5,
  cex.lab = 1.5
)
points(seq(1, 10),
       rep(c_fdr, length(rho)),
       pch = 19 ,cex = 1.5,
       col = 'red')
plot_errorbars(seq(1, 10), rep(c_fdr, length(rho)), rep(c_var, length(rho)))
