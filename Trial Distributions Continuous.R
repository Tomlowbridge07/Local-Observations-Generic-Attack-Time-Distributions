source("General attack solver")

#This file contains pdf and cdf's for testing the general attack solver
#They are all on the support 0 to 3.7

#Increasing

IncreasingLinearTestPDF<-function(x)
{
  if(x>=0 && x<=3.7)
  {
   return(2/(3.7^2)*x)
  }
  else
  {
    return(0)
  }
  
}

IncreasingLinearTestCDF=GenerateCDF(IncreasingLinearTestPDF)

IncreasingQuadraticTestPDF<-function(x)
{
  if(x>=0 && x<=3.7)
  {
   return(3/(3.7^3)*x^2)
  }
  else
  {
    return(0)
  }
}

IncreasingQuadraticTestCDF=GenerateCDF(IncreasingQuadraticTestPDF)

#Decreasing

DecreasingLinearTestPDF<-function(x)
{
  if(x>=0 && x<=3.7)
  {
   return(2/3.7 *(1-x/3.7))
  }
  else
  {
    return(0)
  }
}

DecreasingLinearTestCDF=GenerateCDF(DecreasingLinearTestPDF)

DecreasingQuadraticTestPDF<-function(x)
{
  if(x>=0 && x<=3.7)
  {
   return((5/3.7^3) *(x^2-7.4*x+3.7^2))
  }
  else
  {
   return(0)
  }
}

DecreasingQuadraticTestCDF=GenerateCDF(DecreasingQuadraticTestPDF)

#Uniform

UniformTestPDF<-function(x)
{
  if(x>=0 && x<=3.7)
  {
    return(1/3.7)
  }
  else
  {
    return(0)
  }
}

UniformTestCDF=GenerateCDF(UniformTestPDF)

#Non-monotone, Symmetric type triangular distributions

#Triangular distribution generator (Not intended for use (as just use CDF))
TriangularDistributionPDFGenerator<-function(a,b,c)
{
  TriDisPDF<-function(x)
  {
    if(x<=a)
    {
      return(0)
    }
    else if(x<c)
    {
      return((2*(x-a))/((b-a)*(c-a)))
    }
    else if(x==c)
    {
      return(2/(b-a))
    }
    else if(x<=b)
    {
      return((2*(b-x))/((b-a)*(b-c)))
    }
    else
    {
      return(0)
    }
  }
  return(TriDisPDF)
}

#CDF For Actual use
TriangularDistributionCDFGenerator<-function(a,b,c)
{
  TriDisCDF<-function(x)
  {
    if(x<=a)
    {
      return(0)
    }
    else if(x<=c)
    {
      return(((x-a)^2)/((b-a)*(c-a)))
    }
    else if(x<b)
    {
      return(1-((b-x)^2)/((b-a)*(b-c)))
    }
    else
    {
      return(1)
    }
  }
}


SymmetricTriangularTestCDF=TriangularDistributionCDFGenerator(0,3.7,1.85)

NTypeQuadraticTestPDF<-function(x)
{
  if(x>=0 && x<=3.7)
  {
   return((6/3.7^3) * (3.7*x - x^2)) 
  }
  else
  {
    return(0)
  }
}

NTypeQuadraticTestCDF=GenerateCDF(NTypeQuadraticTestPDF)

RaisedCosineTestPDF<-function(x)
{
  if(x>=0 && x<=3.7)
  {
    return(1/3.7 * (1+cos(pi * ((x-1.85)/1.85)))) 
  }
  else
  {
    return(0)
  }
}

RaisedCosineTestCDF=GenerateCDF(RaisedCosineTestPDF)

#Non-monotone, Un-Symmetric type triangular distributions

LeftSkewTriangularTestCDF=TriangularDistributionCDFGenerator(0,3.7,0.8)

RightSkewTriangularTestCDF=TriangularDistributionCDFGenerator(0,3.7,2.9)

SmallRightSkewTriDisCDF=TriangularDistributionCDFGenerator(0,3.7,2.1)

#ArcSine distribution (BiModal at start and end)

ArcSineTestPDF<-function(x)
{
  if(x>0 && x<3.7)
  {
    return(1/(pi*sqrt((x-0)*(3.7-x))))
  }
  else
  {
    return(0)
  }
}

ArcSineTestCDF=GenerateCDF(ArcSineTestPDF)

TruncatedNormalPDFGenerator<-function(mean,sd,a,b)
{
  TruncNormalPDF<-function(x)
  {
    if(x<a)
    {
      return(0)
    }
    else if(x<=b)
    {
      return((dnorm((x-mean)/sd))/(sd*(pnorm(((b-mean)/sd))-pnorm((a-mean)/sd))))
    }
    else
    {
      return(0)
    }
  }
  return(TruncNormalPDF)
}

TruncatedNormalTest1PDF=TruncatedNormalPDFGenerator(1.7,0.1,0,3.7)
TruncatedNormalTest1CDF=GenerateCDF(TruncatedNormalTest1PDF)
TruncatedNormalTest2PDF=TruncatedNormalPDFGenerator(3.2,0.1,0,3.7)
TruncatedNormalTest2CDF=GenerateCDF(TruncatedNormalTest2PDF)


#Mixing distributions

MixDistributionPDFGenerator<-function(PDF1,PDF2,scaling1)
{
  scaling2=1-scaling1
  MixPDF<-function(x)
  {
    return(scaling1*PDF1(x)+scaling2*PDF2(x))
  }
  return(MixPDF)
}

MixedNormalPDF=MixDistributionPDFGenerator(TruncatedNormalTest1PDF,TruncatedNormalTest2PDF,0.2)
MixedNormalCDF=GenerateCDF(MixedNormalPDF)