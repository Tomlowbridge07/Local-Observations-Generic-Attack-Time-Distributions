source("Truncated Poisson distribution.R")

#This file contains fuctions to solve the problem given g the equilibrium cost

GenerateCDF<-function(AttackTimePDFDistribution,Type="Continuous")
{
  if(Type=="Continuous")
  {
   CDFDistribtion<-function(value)
   {
    return(integrate(Vectorize(AttackTimePDFDistribution),0,value)$value)
   }
   return(CDFDistribtion)
  }
  else if(Type=="Discrete")
  {
    CDFDistribtion<-function(value)
    {
     #Use any value to return pmf matrix
     DiscreteMatrix=AttackTimePDFDistribution(1)
     NumOfSup=ncol(DiscreteMatrix)
     CumProb=vector(length=NumOfSup)
     for(i in 1:NumOfSup)
     {
       if(i==1)
       {
         CumProb[1]=DiscreteMatrix[2,1]
         if(value<DiscreteMatrix[1,1])
         {
           return(0)
         }
         else if(value<DiscreteMatrix[1,2])
         {
           return(CumProb[1])
         }
       }
       else if(i!=1 && i!=NumOfSup)
       {
         CumProb[i]=CumProb[i-1]+DiscreteMatrix[2,i]
         if(value<DiscreteMatrix[1,i])
         {
           return(CumProb[i-1])
         }
         else if(value<DiscreteMatrix[1,i+1])
         {
           return(CumProb[i])
         }
       }
       else
       {
         CumProb[NumOfSup]==CumProb[NumOfSup-1]+DiscreteMatrix[2,NumOfSup]
         if(value<DiscreteMatrix[1,NumOfSup])
         {
           return(CumProb[NumOfSup])
         }
         else if(value>=DiscreteMatrix[1,NumOfSup])
         {
           return(1)
         }
       }
     }
    }
    return(CDFDistribtion)
  }  
}

#Generate the cost to progress matrix
#Note. NumRow is B+1, NumCol=b+1
GenerateSimpleCostToProgressMatrix<-function(NumRow,NumCol,Cost,Lambda,AttackTimeCDFDistribution)
{
  CostToProgressMatrix=matrix(nrow=NumRow,ncol=NumCol)
  for(i in 1:NumRow)
  {
    #We can find the arrivals during this time period
    CostDueToArrivals= Cost * Lambda * integrate(Vectorize(AttackTimeCDFDistribution),i-1,i)$value
    
    for(j in 1:NumCol)
    {
      CostDueToObs= Cost * (j-1) * (AttackTimeCDFDistribution(i)-AttackTimeCDFDistribution(i-1))
      
      CostToProgressMatrix[i,j]=CostDueToArrivals+CostDueToObs
    }
  }
  return(CostToProgressMatrix)
}

#This Generates the cost to progression matrices for arrivals and observed
GenerateCostToProgressMatrix<-function(NumRow,NumCol,Cost,Lambda,AttackTimeCDFDistribution)
{
  CostToProgressMatrix=matrix(nrow=NumRow,ncol=NumCol)
  CostDueToArrivalsMatrix=matrix(nrow=NumRow,ncol=NumCol)
  CostDueToObsMatrix=matrix(nrow=NumRow,ncol=NumCol)
  for(i in 1:NumRow)
  {
    #We can find the arrivals during this time period
    CostDueToArrivals= Cost * Lambda * integrate(Vectorize(AttackTimeCDFDistribution),i-1,i)$value
    
    for(j in 1:NumCol)
    {
      CostDueToObs= Cost * (j-1) * (AttackTimeCDFDistribution(i)-AttackTimeCDFDistribution(i-1))
      
      CostDueToArrivalsMatrix[i,j]=CostDueToArrivals
      CostDueToObsMatrix[i,j]=CostDueToObs
      CostToProgressMatrix[i,j]=CostDueToArrivals+CostDueToObs
    }
  }
  return(list(CostToProgressMatrix=CostToProgressMatrix,CostDueToArrivalsMatrix=CostDueToArrivalsMatrix,CostDueToObsMatrix=CostDueToObsMatrix))
}

#This function returns how long to wait to renew
FindRenewInRow<-function(CurrentRowNum,EquilibriumCost,NextRowRenewIn,CostsToProgress)
{
  CurrentRowRenewIn=vector(length=length(NextRowRenewIn))
 for(j in 1:length(NextRowRenewIn))
 {
   #For element in our row we calculate and check the inequaility (number of g's vs cost)
  CostIfRenew=(NextRowRenewIn[j]+1)*EquilibriumCost
  CostIfNotRenew=0
  for(i in 1:(NextRowRenewIn[j]+1))
  {
   CostIfNotRenew=CostIfNotRenew+CostsToProgress[CurrentRowNum+(i-1),j]
  }
  
  #print(CostIfRenew)
  #print(CostIfNotRenew)
  
  #Compare and decide
  if(CostIfRenew <= CostIfNotRenew)
  {
    #we should renew, set currentRowRenewIn=0
    #print("Chosen to renew")
    CurrentRowRenewIn[j]=0
  }
  else
  {
    #We should wait
    #print("Chosen to wait")
    CurrentRowRenewIn[j]=NextRowRenewIn[j]+1
  }
 }
  return(CurrentRowRenewIn)
}

#This function returns a matrix for each state of how long to wait till renew
FindRenewInMatrix<-function(EquilibriumCost,CostsToProgress)
{
  RenewInMatrix=matrix(nrow=nrow(CostsToProgress),ncol=ncol(CostsToProgress))
  for(rownum in nrow(CostsToProgress):1)
  {
    #print(rownum)
    if(rownum==nrow(CostsToProgress))
    {
      #Preset these to be renew immediately (as we assume g<= c lambda)
      RenewInMatrix[rownum,]=rep(0,ncol(RenewInMatrix))
    }
    else
    {
      #We run the code to figure out the next row up
      RenewInMatrix[rownum,]=FindRenewInRow(rownum,EquilibriumCost,RenewInMatrix[rownum+1,],CostsToProgress)
      #print(RenewInMatrix)
    }
  }
  return(RenewInMatrix)
}

#This function takes a current plan an produces the value of g it would generate
FindEquilibriumValue<-function(RenewInMatrix,CostsToProgress,Lambda,Omega)
{
  #Extract first row of Renew in matrix
  RenewIn=RenewInMatrix[1,]
  CostForCol=vector(length=length(RenewIn))
  ExpectedContributionForCol=vector(length=length(RenewIn))
  ExpectedRenewalLengthContribution=vector(length=length(RenewIn))
  for(j in 1:length(RenewIn))
  {
    #calculate the waiting costs
    if(RenewIn[j]!=0)
    {
      CostForCol[j]=Omega
      for(i in 1:RenewIn[j])
      {
        CostForCol[j]=CostForCol[j]+CostsToProgress[i,j]
      }
    }
    else
    {
      CostForCol[j]=Omega
    }
    
    ExpectedContributionForCol[j]=TruncPoissonPMF(Lambda,length(RenewIn)-1,j) * CostForCol[j]
    ExpectedRenewalLengthContribution[j]=TruncPoissonPMF(Lambda,length(RenewIn)-1,j) * (RenewIn[j]+1)
  }
  ExpectedCostPerRenewal=sum(ExpectedContributionForCol)
  ExpectedRenewalLength=sum(ExpectedRenewalLengthContribution)
  EquilibriumCost=ExpectedCostPerRenewal/ExpectedRenewalLength
  return(EquilibriumCost)
}

#This function generates the equilibrium function as function of omega i.e g(w)
GenerateEquilbriumValue<-function(RenewInMatrix,CostsToProgress,Lambda)
{
  g<-function(Omega)
  {
    return(FindEquilibriumValue(RenewInMatrix,CostsToProgress,Lambda,Omega))
  }
  return(g)
}