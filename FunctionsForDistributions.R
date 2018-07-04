#Example distributions for use

#Functions should return a list of items (pdf value,support type ("d"-discrete,"c"-continuous),support(d-vector of all values,c-vector with min and max values))

TriangularDistributionCreater<-function(Min,Max,Spike)
{
 TriDisFunc<-function(x)
 {
   pdfvalue=NaN
   cdfvalue=NaN
   if(x<Min)
   {
     pdfvalue=0
     cdfvalue=0
   }
   else if(Min<= x && x<= Spike)
   {
     pdfvalue=(2*(x-Min))/((Max-Min)*(Spike-Min))
     cdfvalue=((x-Min)^2)/((Max-Min)*(Spike-Min))
   }
   else if(Spike < x && x<= Max)
   {
     pdfvalue=(2*(Max-x))/((Max-Min)*(Max-Spike))
     cdfvalue=1-((Max-x)^2)/((Max-Min)*(Max-Spike))
   }
   else if(x > Max)
   {
     pdfvalue=0
     cdfvalue=1
   }
   return(list(PDFValue=pdfvalue,CDFValue=,Type="c",Support=c(Min,Max)))
 }
 return(TriDisFunc)
}

UniformDistributionCreater<-function(Min,Max)
{
  UniDisFunc<-function(x)
  {
    pdfvalue=NaN
    cdfvalue=NaN
    if(x<Min)
    {
      pdfvalue=0
      cdfvalue=0
    }
    else if(Min<= x && x<= Max)
    {
      pdfvalue=(1/(Max-Min))
      cdfvalue=(x/(Max-Min))
    }
    else if( x > Max)
    {
      pdfvalue=0
      cdfvalue=1
    }
    return(list(PDFValue=pdfvalue,CDFValue=cdfvalue,Type="c",Support=c(Min,Max)))
  }
  return(UniDisFunc)
}

DiscreteDistributionCreater<-function(Value)
{
  DiscreteDisFunc<-function(x)
  {
    pdfvalue=NaN
    cdfvalue=NaN
    if(x<Value)
    {
      pdfvalue=0
      cdfvalue=0
    }
    else if(x==Value)
    {
      pdfvalue=(1/(Max-Min))
      cdfvalue=1
    }
    else if(x>Value)
    {
      pdfvalue=0
      cdfvalue=1
    }
    return(list(PDFValue=pdfvalue,CDFValue=cdfvalue,Type="d",Support=c(Value)))
  }
  return(DiscreteDisFunc)
}

FunctionPDFRetriver<-function(Function)
{
  PDFValueFunction<-function(x)
  {
    return(Function(x)$PDFValue)
  }
  PDFValueFunction=Vectorize(PDFValueFunction)
  return(PDFValueFunction)
}

#For use with a function with an analytic CDF
FunctionCDFRetriver<-function(Function)
{
  CDFValueFunction<-function(x)
  {
    return(Function(x)$CDFValue)
  }
  CDFValueFunction=Vectorize(CDFValueFunction)
  return(CDFValueFunction)
}