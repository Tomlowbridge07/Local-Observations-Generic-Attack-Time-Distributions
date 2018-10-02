source("General attack solver.R")
source("Iteration Solver.R")

#This file has functions which aim to find an index for the heuristic to follow.
#Note. This heuristic is an approximate only as it uses time steps.


#This function takes in two plans and returns the swapping point omega, by use of the bisection method
FindSwappingPoint<-function(RenewInMatrix1,RenewInMatrix2,CostsToProgress,Lambda,Lowerboundary,Upperboundary,MinErrorTolerance,MaxSteps)
{
  #Form the functions of g
  g1=GenerateEquilbriumValue(RenewInMatrix1,CostsToProgress,Lambda)
  g2=GenerateEquilbriumValue(RenewInMatrix2,CostsToProgress,Lambda)
  f<-function(omega)
  {
    return(g1(omega)-g2(omega))
  }
  
  if(f(Lowerboundary)==0)
  {
    return(Lowerboundary)
  }
  if(f(Upperboundary)==0)
  {
    return(Upperboundary)
  }
  
  #We need to check the boundaries
  #print(f(Lowerboundary))
  #print(f(Upperboundary))
  stopifnot(sign(f(Lowerboundary))!=sign(f(Upperboundary)))
  stopifnot(f(Upperboundary)>f(Lowerboundary))
  
  #Repeat the bisection process until tolerance or number of steps is reached
  tolerance=MinErrorTolerance+1
  Steps=0
  a=Lowerboundary
  b=Upperboundary
  while(Steps<MaxSteps && tolerance>MinErrorTolerance)
  {
    midpoint=(a+b)/2
    
    #evaluate mid point
    midpointvalue=f(midpoint)
    
    #if it is zero return
    if(f(midpoint)==0)
    {
      return(midpoint)
    }
    
    #see what to replace
    if(sign(f(midpoint))==sign(f(a)))
    {
      a<-midpoint
    }
    else
    {
      b<-midpoint
    }
    
    #update steps and tolerance
    Steps=Steps+1
    tolerance=(b-a)/2
    
  }
  return(midpoint)
}

#Return w* for every state. We choose to ignore later ones.
FindIndexMatrix<-function(OmegaStepSize,AttackTimeDistribution,B,b,Cost,Lambda,MinTolerance,MaxSteps,TypeOfAttackTimeDis="CDF")
{
  #We'll generate the cost to progress matrix for later use
  CostToProgress=GenerateCostToProgressMatrix(B+1,b+1,Cost,Lambda,AttackTimeDistribution)
  
  #We first solve to the boundary
  print("Solving for each omega")
  SolvedOE=SolveForMultipleOmegaUntilBoundary(OmegaStepSize,AttackTimeDistribution,B,b,Cost,Lambda,TypeOfAttackTimeDis="CDF")
  
  #Retrive the plans
  Plans=SolvedOE$Plans
  PlansChangingPoints=SolvedOE$PlansChangingPoints
  
  print("Creating Index Matrix")
  #Index matrix for w*'s to be put into.
  IndexMatrix=matrix(nrow=B+1,ncol=b+1)
  
  #Then for each plan change we will find what states change and what.
  #Note. The last Plan is past the boundary
  for(i in 1:(length(Plans)-1))
  {
    #Find w*
    OmegaStar=FindSwappingPoint(Plans[[i]],Plans[[i+1]],CostToProgress,Lambda,PlansChangingPoints[i],PlansChangingPoints[i+1],MinTolerance,MaxSteps)
    
    #print(OmegaStar)
    
    #Find what elements have changed in the plan
    for(element in 1:length(Plans[[i]]))
    {
      if(Plans[[i]][element]==0 && Plans[[i+1]][element]!=0)
      {
        IndexMatrix[element]=OmegaStar
      }
    }
  }
  
  #We pass the boundary hit value on
  BoundaryHitValue=SolvedOE$BoundaryHitValue
  
  return(list(IndexMatrix=IndexMatrix,BoundaryHitValue=BoundaryHitValue))
}

#This function alters the index matrix to have a 
AlterIndexMatrix<-function(IndexMatrix,IndexCap,IgnoreIndexibility=F)
{
  AlteredIndexMatrix=IndexMatrix
  if(IgnoreIndexibility==F)
  {
   #Not Ignoring indexibility we must remove all values which are less on the way down.
   #We must also remove any values set below
    
    
    #Each column when it hits NA should be NA for the rest of the column.
    for(j in 1:ncol(IndexMatrix))
    {
      #go down the column to the first NA, then replace all others with NA
      IsNAMode=F
      for(i in 1:nrow(IndexMatrix))
      {
        if(is.na(IndexMatrix[i,j]))
        {
          IsNAMode=T
        }
        if(IsNAMode==T)
        {
          AlteredIndexMatrix[i,j]=NA
        }
      }
    }
    
     
  }
  else
  {
    #Ignoring indexibility we will keep all current indices. So we have no additional work
  }
  
  #then we replace all NA with our cap value (boundary)
  for(i in 1:nrow(AlteredIndexMatrix))
  {
    for(j in 1:ncol(AlteredIndexMatrix))
    {
      if(is.na(AlteredIndexMatrix[i,j]))
      {
        AlteredIndexMatrix[i,j]=IndexCap
      }
      
    }
  }
  
  return(AlteredIndexMatrix)
}


