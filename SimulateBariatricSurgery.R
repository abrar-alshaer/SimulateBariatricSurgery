#Abrar Al-Shaer, Raz Shaikh Lab
#6/28/19
#This program was adapted from Dr. Anthony Fodor, UNCC Bioinformatics Dept.

# this power simulation assumes a normal distribution 
# and then applies a t-test.

rm(list=ls())
sampleSizeForEachGroup <- 50
numHypotheses <- 31 #number of measured markers
numSimulations <- 1000
fractionTruePositives <- 0.25 #assuming 25% of the markers can be significant 
effectSize <- 0.75 #modeling slightly above moderate effect size

pValues <- vector()
truePositives <- vector()
power_bonf_ten <- vector()
power_bonf_five <- vector()
power_bh_ten <- vector()
power_bh_five <- vector()

for (k in 1:numSimulations) {
  numUnder <- 0
  numBonferroniUnder <-0
  numBonf <- 0
  numTruePositves <-0
  bhUnder <- 0
  benjHochUnder <- 0
  for( j in 1 : numHypotheses ) 
  {
    isATruePositive <- ( runif(1) <= fractionTruePositives )
    
    if( isATruePositive ) 	
      numTruePositves= numTruePositves + 1
    
    data<-vector()
    data2<-vector()
    
    for( i in 1 : sampleSizeForEachGroup)
    {
      data[i] = rnorm(1) #add regardless if it's a true positive
      
      if( isATruePositive  ) #add to data2 ONLY if it is a true positive
      {
        data2[i] = rnorm(1, mean=effectSize)
      }
      else
      {
        data2[i] = rnorm(1)
      }
    }
    
    pValues[j] <- t.test(data,data2)$p.value #is there is a significant difference between data and data2
    truePositives[j] <- isATruePositive
    
    if( isATruePositive  & pValues[j] < 0.05 )
      numUnder = numUnder + 1
    
    # if( isATruePositive  & pValues[j] < 0.05 / numHypotheses )
    #   numBonferroniUnder = numBonferroniUnder + 1	
    
  } 
  
  numTruePositves
  numUnder/numTruePositves
  
  myFrame <- data.frame(pValues,truePositives  )
  myFrame<- myFrame[order(myFrame$pValues),]
  myFrame$pAdjust <- p.adjust(myFrame$pValues, method="BH")
  myFrame$Bonf <- p.adjust(myFrame$pValues, method="bonferroni")
  
  for( a in 1:nrow(myFrame))
  {
    if( myFrame$truePositives[a] & myFrame$pAdjust[a] < 0.10  )
      bhUnder <- bhUnder + 1
    
    if( myFrame$truePositives[a] & myFrame$pAdjust[a] < 0.05  )
      benjHochUnder <- benjHochUnder + 1
    
    if( myFrame$truePositives[a] & myFrame$Bonf[a] < 0.10  )
      numBonferroniUnder = numBonferroniUnder + 1
    
    if( myFrame$truePositives[a] & myFrame$Bonf[a] < 0.05  )
      numBonf = numBonf + 1
  } 
  
  power_bonf_ten[k] <- numBonferroniUnder/numTruePositves
  power_bonf_five[k] <- numBonf/numTruePositves
  #print(numBonferroniUnder/numTruePositves)
  
  power_bh_ten[k] <- bhUnder / numTruePositves
  power_bh_five[k] <- benjHochUnder / numTruePositves
  #print(bhUnder / numTruePositves) 
}

mean(na.omit(power_bonf_ten))
mean(na.omit(power_bonf_five))
mean(na.omit(power_bh_ten))
mean(na.omit(power_bh_five))

