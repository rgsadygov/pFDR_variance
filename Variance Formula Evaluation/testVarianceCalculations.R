library("ggplot2")
source("./pFDR_and_Variance.r")
source("./cleaned_pFDR_and_Var.r")


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
  
  columnNames <- c("alpha", "beta", "m0", "m1", "FDR_Integration", "VarFDP_Binomal", "VarFDP_Integration", "relativeDifference")
  toleranceFrame <- data.frame(matrix(ncol=length(columnNames), nrow=0))
  colnames(toleranceFrame) <- columnNames #generate dataframe to store parameters and calculate values
  
  for(alpha in alpha_range) {
    for(beta in beta_range) {
      for(m0 in m0_range) {
        for(m1 in m1_range) {
          
          VarQSummation <- calculateVarQSummation(alpha, beta, m0, m1) #Variance of FDP from binomial approach
          VarQIntegral <- pFDR_variance(alpha, beta, m0, m1) #Variance of FDP from integration approach
          relTolerance <- (VarQSummation - VarQIntegral)/VarQSummation #determine relative difference with respect to binomial approach
          
          FDRintegral <- Integral_pFDR(alpha, beta, m0, m1) #expectation of FDP aka FDR
          newRow = c(alpha, beta, m0, m1, FDRintegral, VarQSummation, VarQIntegral, relTolerance)
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

