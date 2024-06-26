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



#=================================================================
#===================FDP Integral approach=========================
#=================================================================
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

#calculates the pFDR using the integration approach
pFDR_IntegralApproach <- function(alpha, beta, m0, m1, mid_step= -Inf)
{
  N_grid    = 10;
  mid_step = 1;
  N_crit = m0 + m1;
  while (mid_step == 1)
  {
    mid_step = critical_point_1d(alpha, beta, m0, m1, N_crit);
    N_crit = 10 * N_crit;
  }
  
  if(mid_step > 10)
  {
    mid_step = mid_step - 10;
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
  
  return (pFDR_theoretical)
}
 
#=================================================================
#===================END FDP Integral approach=========================
#=================================================================


#
# The function, compute_fdp) computes the false discovery proportion given, rho(correlation),
# alpha (significane level of the null), and beta (the power of the test) 
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
  
  print(c(
    'rho ',
    rho,
    'FDP ',
    fdp,
    'Integral FDP ',
    pFDR_IntegralApproach(alpha, beta, n.sim, n.sim)
  ))
  return(fdp)
}


#======================================================================

# run the simulation for diffrent Correlation values
set.seed(10)
rho = c(-0.99,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.2, 0.1,-0.0001)
alpha = 0.1
beta = 0.9
n.sim <- 10000
fdps = c(rep(-1, length(rho)))
for (i in seq(1, length(rho))) {
  fdps[i] = compute_fdp(rho[i], alpha, beta, n.sim)
}

plot(
  rho,
  fdps,
  pch = 19,
  ylim = c(min(fdps) - 0.05, max(fdps) + 0.05),
  cex = 1.5,
  cex.axis = 1.5,
  cex.lab = 1.5,
  ylab = 'FDP',
  xlab = 'Correlation coefficient'
)
points(
  rho,
  rep(pFDR_IntegralApproach(alpha, beta, n.sim, n.sim), length(rho)),
  col = 'red' ,
  cex = 1.5,
  pch = 19
)

legend(
  "topright",
  legend = c('FDP', 'FDR from integral formula'),
  col = c('black', 'red'),
  pch = 19,
  bty = "n"
)
