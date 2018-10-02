source("Index Implementation.R")

DeterministicCostEvaluationOfPath<-function(Path,n,sVec,vVec,CostVec,LambdaVec,BVec,bVec,CostToProgressList)
{
  #We follow the path for as many steps as it has
  i=1
  SumOfCost=0
  while(Path[i]!=0 && i<=length(Path))
  {
    #Add actions cost
    SumOfCost=SumOfCost+CostOfAction(c(sVec,vVec),Path[i],n,CostToProgressList)
    
    #Evolve DETERMINISTICALLY
    sVec=NewSState(sVec,Path[i],BVec)
    vVec=NewMeanVState(vVec,sVec,Path[i],BVec,bVec,LambdaVec)

    i=i+1
  }
  
  return(list(NewMeanSVState=c(sVec,vVec),Average=SumOfCost/(i-1),Overall=SumOfCost))
}

#This function finds the cost while we neglect to visit, until all states are full
DecayToEndCost<-function(n,sVec,vVec,CostVec,LambdaVec,bVec,BVec,AttackCDFVec,CostToProgressList=NULL)
{
  if(is.null(CostToProgressList))
  {
    CostToProgressList=CreateCostToProgressList(BVec,bVec,CostVec,LambdaVec,AttackCDFVec)
  }
  
  TransitionaryCosts=vector(length=0)
  #While we haven't reached the end we sum some costs
  while(!all(sVec==(BVec+1)))
  {
    NeglectingCost=CostToNeglect(c(sVec,vVec),n,CostToProgressList)
    NewState=NewSVStateIfNeglect(c(sVec,vVec),BVec,bVec,LambdaVec)
    sVec=NewState[1:n]
    vVec=NewState[(n+1):(2*n)]
    TransitionaryCosts=c(TransitionaryCosts,NeglectingCost)
  }
  
  if(length(TransitionaryCosts)==0)
  {
    return(0)
  }
  else
  {
   Sum=sum(TransitionaryCosts)
   Avg=Sum/length(TransitionaryCosts)
   return(Avg)
  }
}

#Here we sum the indices up over the number of steps for all paths of the Number of steps
MultiStepBenefitHeuristic<-function(NoSteps,n,AdjacencyMatrix,IndexList,sVec,vVec,CostVec,LambdaVec,BVec,bVec,AttackCDFVec,CostToProgressList=NULL,PrintOutput=FALSE)
{
  #Generate cost to progress list if not given
  if(is.null(CostToProgressList))
  {
    CostToProgressList=CreateCostToProgressList(BVec,bVec,CostVec,LambdaVec,AttackCDFVec)
  }
  
  if(PrintOutput)
  {
   print("Heuristic is Run")
  }
  
  
  Paths=matrix(0,nrow=0,ncol=NoSteps)
  #Note. We will be using mean v Evolution for use with the index
  EvolvedStates=matrix(0,nrow=0,ncol=2*n)
  BenefitForPath=vector(length=0)
  #Note. Best Path for Step may use multiple for steps if there are mutliple option for that length of path
  BestPathforStep=matrix(0,nrow=0,ncol=NoSteps)
  
  for(Step in 1:NoSteps)
  {
    if(Step==1)
    {
      #We run the initial set up
      CurrentNode=which.min(sVec)
      
      #Actions
      Actions=AdjacencyMatrix[CurrentNode,]
      
      BenefitForAction=vector(length=length(Actions))
      for(action in 1:length(BenefitForAction))
      {
        BenefitForAction[action]=0
        if(Actions[action]==1)
        {
          BenefitForAction[action]=IndexFromListReader(action,sVec[action],vVec[action],IndexList)
        }
      }
      
      #Form Paths
      for(action in 1:n)
      {
        if(Actions[action]==1)
        {
          Paths=rbind(Paths,c(action,rep(0,NoSteps-1)))
          NewsVec=NewSState(sVec,action,BVec)
          EvolvedStates=rbind(EvolvedStates,c(NewsVec,NewMeanVState(vVec,NewsVec,action,BVec,bVec,LambdaVec)))
          print(EvolvedStates)
          BenefitForPath=c(BenefitForPath,BenefitForAction[action])

        }
      }

      if(PrintOutput)
      {
       print(Paths)
       print(BenefitForPath)
      }
      
      #Identify the maximal elements
      MaximalElements=which(BenefitForPath==max(BenefitForPath))

      
      #We now store all those that provide maximal benefit
      for(MaximalElementsNo in 1:length(MaximalElements))
      {
        BestPathforStep=rbind(BestPathforStep,Paths[MaximalElements[MaximalElementsNo],])
      }
      
      
      
      
      # print("I have chosen the path")
      # print(BestPath)
    }
    else
    {
      
      #We need to copy and expand each
      OldPaths=Paths
      OldEvolvedStates=EvolvedStates
      OldBenefitForPath=BenefitForPath
      
      Paths=matrix(0,nrow=0,ncol=NoSteps)
      EvolvedStates=matrix(0,nrow=0,ncol=2*n)
      BenefitForPath=vector(length=0)
      
      #We now look at expanding each
      for(row in 1:nrow(OldPaths))
      {
        #for each row we expand it to allow all possible actions
        Actions=AdjacencyMatrix[OldPaths[row,Step-1],]
        
        BenefitForAction=vector(length=length(Actions))
        for(action in 1:length(BenefitForAction))
        {
          BenefitForAction[action]=0
          if(Actions[action]==1)
          {
            print("Evolved state is")
            print(c(OldEvolvedStates[row,1:n][action],OldEvolvedStates[row,(n+1):(2*n)][action]))
            print(CostToProgressList)
            print(IndexFromListReader(action,OldEvolvedStates[row,1:n][action],OldEvolvedStates[row,(n+1):(2*n)][action],IndexList))
           BenefitForAction[action]=IndexFromListReader(action,OldEvolvedStates[row,1:n][action],OldEvolvedStates[row,(n+1):(2*n)][action],IndexList)
          }
        }

        
        #Form Paths
        for(action in 1:n)
        {
          if(Actions[action]==1)
          {
            Paths=rbind(Paths,c(OldPaths[row,1:(Step-1)],action,rep(0,NoSteps-Step)))
            NewsVec=NewSState(OldEvolvedStates[row,1:n],action,BVec)
            EvolvedStates=rbind(EvolvedStates,c(NewsVec,NewMeanVState(OldEvolvedStates[row,(n+1):(2*n)],NewsVec,action,BVec,bVec,LambdaVec)))
            BenefitForPath=c(BenefitForPath,OldBenefitForPath[row]+BenefitForAction[action])
          }
        }
        
      }
      # print(paste("I am about to compare all paths of length ",toString(Step)))
      if(PrintOutput)
      {
        print(Paths)
        print(BenefitForPath)
      }
      #Identify the maximal elements
      MaximalElements=which(BenefitForPath==max(BenefitForPath))
      # print("Printing maximal elements")
      # print(MaximalElements)

      
      #We now store all those that provide maximal benefit
      for(MaximalElementsNo in 1:length(MaximalElements))
      {
        BestPathforStep=rbind(BestPathforStep,Paths[MaximalElements[MaximalElementsNo],])
      }
      
      # print("I have chosen the path")
      # print(BestPath)

    }
    
  }
  
  #For each look ahead step we have a path
  #print(BestPathforStep)
  NumberOfPathsToConsider=nrow(BestPathforStep)
  AverageCostForPath=vector(length=NumberOfPathsToConsider)
  #We now need to see how good they perform
  for(i in 1:NumberOfPathsToConsider)
  {
    #We compute the average cost of following such a strategy to decide which paths to pick
    #We use determinsitic evolution to the mean state in v
    DeterministicGen=DeterministicCostEvaluationOfPath(BestPathforStep[i,],n,sVec,vVec,CostVec,LambdaVec,BVec,bVec,CostToProgressList)
    AverageCostForPath[i]=DeterministicGen$Average
    #For each Path we have have end point
    EndOfPathState=DeterministicGen$NewMeanSVState
    #From this work out the average cost to decay
    AverageCostForPath[i]=AverageCostForPath[i]+DecayToEndCost(n,EndOfPathState[1:n],EndOfPathState[(n+1):(2*n)],CostVec,LambdaVec,bVec,BVec,AttackCDFVec,CostToProgressList)
  }
  # print("about to print paths and determinisitic cost of paths")
  if(PrintOutput)
  {
   print(BestPathforStep)
   print(AverageCostForPath)
  }

  #Identify the maximal elements
  MinimalElements=which(AverageCostForPath==min(AverageCostForPath))
  #We now choose one at random
  #ChosenMin=MinimalElements[sample(1:length(MinimalElements),1)]
  #We will just pick the minimum to avoid confusion
  ChosenMin=MinimalElements[1]
  OverallBestPath=BestPathforStep[ChosenMin,]
  return(OverallBestPath[1])
}


#Here we sum the indices up over the number of steps for all paths of the Number of steps
MultiStepPenaltyHeuristic<-function(NoSteps,n,AdjacencyMatrix,IndexForNodeFunction,sVec,vVec,CostVec,LambdaVec,bVec,xVec,vMaxVec=NULL,PrintOutput=FALSE)
{
  if(PrintOutput)
  {
    print("Heuristic is Run")
  }
  if(is.null(vMaxVec))
  {
    #Create vMaxVec
    vMaxVec=CreateVMaxVector(n,LambdaVec,bVec,xVec)
  }
  
  BVec=ceiling(xVec)
  
  
  Paths=matrix(0,nrow=0,ncol=NoSteps)
  #Note. We will be using mean v Evolution for use with the index
  EvolvedStates=matrix(0,nrow=0,ncol=2*n)
  PenaltyForPath=vector(length=0)
  BestPathforStep=matrix(0,nrow=0,ncol=NoSteps)
  
  for(Step in 1:NoSteps)
  {
    if(Step==1)
    {
      #We run the initial set up
      CurrentNode=which.min(sVec)
      
      #Actions
      Actions=AdjacencyMatrix[CurrentNode,]
      
      NodeIndexes=IndicesForNodes(n,IndexForNodeFunction,sVec,vVec,CostVec,LambdaVec,bVec,xVec,vMaxVec)
      
      IndexForActions= Actions * NodeIndexes
      
      PenaltyForAction=rep(sum(IndexForActions),length(Actions))-IndexForActions
      
      #Form Paths
      for(action in 1:n)
      {
        if(Actions[action]==1)
        {
          Paths=rbind(Paths,c(action,rep(0,NoSteps-1)))
          NewsVec=NewSState(sVec,action,BVec)
          EvolvedStates=rbind(EvolvedStates,c(NewsVec,NewMeanVState(vVec,NewsVec,action,BVec,bVec,LambdaVec)))
          PenaltyForPath=c(PenaltyForPath,PenaltyForAction[action])
        }
      }
      if(PrintOutput)
      {
        print(Paths)
        print(PenaltyForPath)
      }
      

      
      #Identify the maximal elements
      MinimalElements=which(PenaltyForPath==min(PenaltyForPath))
      #We now choose one at random
      #ChosenMin=MinimalElements[sample(1:length(MinimalElements),1)]
      
      #We now store all those that provide maximal benefit
      for(MinimalElementsNo in 1:length(MinimalElements))
      {
        BestPathforStep=rbind(BestPathforStep,Paths[MinimalElements[MinimalElementsNo],])
      }
    }
    else
    {
      #We need to copy and expand each
      OldPaths=Paths
      OldEvolvedStates=EvolvedStates
      OldPenaltyForPath=PenaltyForPath
      
      Paths=matrix(0,nrow=0,ncol=NoSteps)
      EvolvedStates=matrix(0,nrow=0,ncol=2*n)
      PenaltyForPath=vector(length=0)
      
      #We now look at expanding each
      for(row in 1:nrow(OldPaths))
      {
        #for each row we expand it to allow all possible actions
        Actions=AdjacencyMatrix[OldPaths[row,Step-1],]
        
        NodeIndexes=IndicesForNodes(n,IndexForNodeFunction,OldEvolvedStates[row,1:n],OldEvolvedStates[row,(n+1):(2*n)],CostVec,LambdaVec,bVec,xVec,vMaxVec)
        
        IndexForActions= Actions * NodeIndexes
        
        PenaltyForAction=rep(sum(IndexForActions),length(Actions))-IndexForActions
        
        #Form Paths
        for(action in 1:n)
        {
          if(Actions[action]==1)
          {
            Paths=rbind(Paths,c(OldPaths[row,1:(Step-1)],action,rep(0,NoSteps-Step)))
            NewsVec=NewSState(OldEvolvedStates[row,1:n],action,BVec)
            EvolvedStates=rbind(EvolvedStates,c(NewsVec,NewMeanVState(OldEvolvedStates[row,(n+1):(2*n)],NewsVec,action,BVec,bVec,LambdaVec)))
            PenaltyForPath=c(PenaltyForPath,OldPenaltyForPath[row]+PenaltyForAction[action])
          }
        }
        
      }
      if(PrintOutput)
      {
        print(Paths)
        print(PenaltyForPath)
      }
      #Identify the maximal elements
      MinimalElements=which(PenaltyForPath==min(PenaltyForPath))
      #We now choose one at random
      #ChosenMin=MinimalElements[sample(1:length(MinimalElements),1)]
      #We now store all those that provide maximal benefit
      
      #We now store all those that provide maximal benefit
      for(MinimalElementsNo in 1:length(MinimalElements))
      {
        BestPathforStep=rbind(BestPathforStep,Paths[MinimalElements[MinimalElementsNo],])
      }
    }
    
  }
  
  #For each look ahead step we have a path
  #print(BestPathforStep)
  AverageCostForPath=vector(length=NoSteps)
  #We now need to see how good they perform
  for(i in 1:NoSteps)
  {
    #We compute the average cost of following such a strategy to decide which paths to pick
    #We use determinsitic evolution to the mean state in v
    DeterministicGen=DeterministicCostEvaluationOfPath(BestPathforStep[i,],n,sVec,vVec,CostVec,LambdaVec,bVec,xVec,vMaxVec,CostToProgressList)
    AverageCostForPath[i]=DeterministicGen$Average
    #For each Path we have have end point
    EndOfPathState=DeterministicGen$NewMeanSVState
    #From this work out the average cost to decay
    AverageCostForPath[i]=AverageCostForPath[i]+DecayToEndCost(n,EndOfPathState[1:n],EndOfPathState[(n+1):(2*n)],CostVec,LambdaVec,bVec,xVec)
    
    
  }
  #print(AverageCostforPath)
  if(PrintOutput)
  {
   print(BestPathforStep)
   print(AverageCostForPath)
  }

  #Identify the maximal elements
  MinimalElements=which(AverageCostForPath==min(AverageCostForPath))
  #We now choose one at random
  #ChosenMin=MinimalElements[sample(1:length(MinimalElements),1)]
  ChosenMin=MinimalElements[1]
  OverallBestPath=BestPathforStep[ChosenMin,]
  return(OverallBestPath[1])
}


StartingNodeHeuristic<-function(n,IndexForNodeFunction,CostVec,LambdaVec,bVec,xVec,vMaxVec=NULL)
{
  if(is.null(vMaxVec))
  {
    #Create vMaxVec
    vMaxVec=CreateVMaxVector(n,LambdaVec,bVec,xVec)
  }
  
  #Elapsed state, we will consider the indices when we have not visited in a very long time.
  StartingSVec=ceiling(xVec)+1
  #and we will assume that all v's are in their mean state.
  StartingVVec=vector(length=n)
  for(i in 1:n)
  {
    # StartingVVec[i]=TruncPoissionMean(LambdaVec[i],bVec[i])
    StartingVVec[i]=0
  }

  #To decide a starting node we work out the index for all nodes (when S is maximum) and pick the biggest
  NodeIndexes=IndicesForNodes(n,IndexForNodeFunction,StartingSVec,StartingVVec,CostVec,LambdaVec,bVec,xVec,vMaxVec)
  print("Here")
  BestStart=which.max(NodeIndexes)
  return(BestStart)
}

HeuristicPolicy<-function(HeuristicDepth,HeuristicFunction,n,AdjacencyMatrix,IndexForNodeFunction,CostVec,LambdaVec,bVec,xVec,StateSpace=NULL,vMaxVec=NULL,PrintOutput=FALSE)
{
  BVec=ceiling(xVec)
  Policy=list()
  if(is.null(StateSpace))
  {
    StateSpace=CreateSVStates(n,BVec,bVec)
  }
  
  #For each state we will apply our algorithm to get a policy
  for(StateNumber in 1:nrow(StateSpace))
  {
    #For each state we will find what the Heurisitic tells us to do
    State=StateSpace[StateNumber,]
    if(PrintOutput)
    {
      print(State)
    }
    MoveToNode=HeuristicFunction(HeuristicDepth,n,AdjacencyMatrix,IndexForNodeFunction,State[1:n],State[(n+1):(2*n)],CostVec,LambdaVec,bVec,xVec,vMaxVec)
    
    #We record in a list the policy
    Policy[[StateNumber]]=MoveToNode
  }
  return(Policy)
}