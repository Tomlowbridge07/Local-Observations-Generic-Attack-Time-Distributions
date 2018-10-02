source("General attack solver.R")
library(ggplot2)
library(reshape)
library(tikzDevice)

IterationSolver<-function(Omega,AttackTimeDistribution,B,b,Cost,Lambda,TypeOfAttackTimeDis="CDF")
{
 if(TypeOfAttackTimeDis=="PDF")
 {
   #Conversion is required
   AttackTimeDistribution=GenerateCDF(AttackTimePDFDistribution)
 }
  
 #Generate cost matrix
 CostToProgressMatrix=GenerateCostToProgressMatrix(B+1,b+1,Cost,Lambda,AttackTimeDistribution)
 
 #Then we peform iterations to find g with an initial guess of c * lambda
 InitialGuess= Cost * Lambda
 EvaluatedEquilibrium=InitialGuess
 MinTolerance=0.0001
 MaxSteps=1000
 Step=1
 Tolerance=MinTolerance+1
 CurrentPlan=FindRenewInMatrix(InitialGuess,CostToProgressMatrix)
 while(Tolerance>MinTolerance && Step<MaxSteps)
 {
   OldEquilibriumValue=EvaluatedEquilibrium
   EvaluatedEquilibrium=FindEquilibriumValue(CurrentPlan,CostToProgressMatrix,Lambda,Omega)
   CurrentPlan=FindRenewInMatrix(EvaluatedEquilibrium,CostToProgressMatrix)
   
   #Now see if equilibrium value has changed (Equally we could see if plan has changed ??)
   
   Tolerance=abs(EvaluatedEquilibrium-OldEquilibriumValue)
   Step=Step+1
 }
 return(list(Plan=CurrentPlan,EquilibriumValue=EvaluatedEquilibrium))
}

SolveForMultipleOmega<-function(OmegaMin,OmegaMax,OmegaSteps,AttackTimeDistribution,B,b,Cost,Lambda,TypeOfAttackTimeDis="CDF")
{
 #For each omega we run the code, store g for that omega
 OmegaIncrease=(OmegaMax-OmegaMin)/OmegaSteps
 CurrentPlan=matrix(rep(-1,(B+1)*(b+1)),nrow=B+1,ncol=b+1)
 OmegaEquilibrium=matrix(nrow=2,ncol=OmegaSteps+1)
 FullPlans=list(length=OmegaSteps+1)
 PlanChanging=vector(length=0)
 PlanChangingEqValues=vector(length=0)
 Plans=list()
 PlanChangingCounter=1
 BoundaryHit=F
 
 for(i in 1:(OmegaSteps+1))
 {
   Omega=OmegaMin+(i-1)*OmegaIncrease
   #Run iteration solver to find plan and g
   Solved=IterationSolver(Omega,AttackTimeDistribution,B,b,Cost,Lambda,TypeOfAttackTimeDis)
   g=Solved$EquilibriumValue
   Plan=Solved$Plan
   FullPlans[[i]]=Plan
   OmegaEquilibrium[1,i]=Omega
   OmegaEquilibrium[2,i]=g
   
   #See if the plan changes
   if(!all(CurrentPlan==Plan))
   {
     PlanChanging=c(PlanChanging,Omega)
     Plans[[PlanChangingCounter]]=Plan
     PlanChangingEqValues=c(PlanChangingEqValues,g)
     PlanChangingCounter=PlanChangingCounter+1
   }
   CurrentPlan=Plan
   
   #Find when it hits the boundary
   if(BoundaryHit==F && g> Cost * Lambda)
   {
     BoundaryHitValue=Omega
     BoundaryHit=T
   }
 }
 
 #If the boundary hit is still false then we return NUll
 if(BoundaryHit==F)
 {
   BoundaryHitValue=NULL
 }
 
 return(list(OmegaEquilibriumMatrix=OmegaEquilibrium,FullPlans=FullPlans,PlansChangingPoints=PlanChanging,
             PlanChangingEqValues=PlanChangingEqValues,Plans=Plans,BoundaryHitValue=BoundaryHitValue))
}

#Similiar to above but runs until the boundary is hit
SolveForMultipleOmegaUntilBoundary<-function(OmegaStepSize,AttackTimeDistribution,B,b,Cost,Lambda,TypeOfAttackTimeDis="CDF")
{
  #For each omega we run the code, store g for that omega
  OmegaIncrease=OmegaStepSize
  CurrentPlan=matrix(rep(-1,(B+1)*(b+1)),nrow=B+1,ncol=b+1)
  OmegaEquilibrium=matrix(nrow=2,ncol=0)
  FullPlans=list()
  PlanChanging=vector(length=0)
  PlanChangingEqValues=vector(length=0)
  Plans=list()
  PlanChangingCounter=1
  BoundaryHit=F
  Omega=0
  i=1
  
  while(BoundaryHit==F)
  {
    #Run iteration solver to find plan and g
    Solved=IterationSolver(Omega,AttackTimeDistribution,B,b,Cost,Lambda,TypeOfAttackTimeDis)
    g=Solved$EquilibriumValue
    Plan=Solved$Plan
    FullPlans[[i]]=Plan
    i=i+1
    OmegaEquilibrium=cbind(OmegaEquilibrium,matrix(c(Omega,g),nrow=2,ncol=1))
    
    #See if the plan changes
    if(!all(CurrentPlan==Plan))
    {
      PlanChanging=c(PlanChanging,Omega)
      Plans[[PlanChangingCounter]]=Plan
      PlanChangingEqValues=c(PlanChangingEqValues,g)
      PlanChangingCounter=PlanChangingCounter+1
    }
    CurrentPlan=Plan
    
    #Find when it hits the boundary
    if(BoundaryHit==F && g> Cost * Lambda)
    {
      BoundaryHitValue=Omega
      BoundaryHit=T
    }
    
    Omega=Omega+OmegaIncrease
  }
  
  return(list(OmegaEquilibriumMatrix=OmegaEquilibrium,FullPlans=FullPlans,PlansChangingPoints=PlanChanging,
              PlanChangingEqValues=PlanChangingEqValues,Plans=Plans,BoundaryHitValue=BoundaryHitValue)) 
}

PlotOmegaEquilibrium<-function(OmegaEquilibriumMatrix,PlanChanging,Cost,Lambda)
{
  XCoordinates=OmegaEquilibriumMatrix[1,]
  YCoordinates=OmegaEquilibriumMatrix[2,]
  YBoundary=rep(Cost*Lambda,length(XCoordinates))
  XYCoordinatedataframe=data.frame(XCoordinates)
  XYCoordinatedataframe<-cbind(XYCoordinatedataframe,YCoordinates)
  XYCoordinatedataframe<-cbind(XYCoordinatedataframe,YBoundary)
  print(XYCoordinatedataframe)
  DataFrame<-XYCoordinatedataframe
  print(DataFrame)
  MeltedDataFrame<-melt(DataFrame,id="XCoordinates")
  print(MeltedDataFrame)
  
  XPointChanging=PlanChanging
  print(XPointChanging)
  YPointChanging=rep(0,length(XPointChanging))
  print(YPointChanging)
  for(i in 1:length(XPointChanging))
  {
   YPointChanging[i]=YCoordinates[which(XCoordinates==XPointChanging[i])]
  }
  print(YPointChanging)
  XYPointdataframe<-data.frame(XPointChanging)
  XYPointdataframe<-cbind(XPointChanging,YPointChanging)
  MeltedDataFrame2<-melt(as.data.frame(XYPointdataframe),id="XPointChanging")
  print(XYPointdataframe)
  print(MeltedDataFrame2)
  
  Plot<-ggplot(MeltedDataFrame,aes(x=XCoordinates,y=value,color=variable),show.legend='True')+geom_line()+
    geom_point(aes(x=XPointChanging,y=YPointChanging),data=MeltedDataFrame2,inherit.aes = F)+
    
  print(Plot)
}