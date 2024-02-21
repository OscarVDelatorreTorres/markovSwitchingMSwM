
MSwMAnalysis2K=function(Data){
  
  TstartTime=Sys.time()
  mod=lm("Return~1",data=Data)
  
  #cat("\f")
  print("Estimating 2-regime MS model...")
  startTime=Sys.time()
  mod.mswm=msmFit(mod,k=2,p=0,sw=c(T,T),control=list(parallel=F))
  endTime=Sys.time()
  estElaptTime=round((endTime-startTime)/60,4)
  
  regimenId=rank(mod.mswm@std)

  cat("\f")
  print("Creating output data...")    
  # Records the observed return:
  
  tempDataFrame=data.frame(Date=as.character(tail(Data$Date,1)),
                           Value=tail(Data$Return,1),
                           Ticker="observedReturn")
  
  outPutDataFrame=tempDataFrame
  
  # Location and scale parameters:
  # Location:
  MSparameter=matrix(0,2,1)
  
  for (a in 1:2){
    MSparameter[a]=mod.mswm@Coef$`(Intercept)`[regimenId[a]]
  }
  
  tempDataFrame=data.frame(
    Date=matrix(as.character(tail(Data$Date,1),4,1)),
    Value=MSparameter,
    Ticker=c("expReturnReg1","expReturnReg2")
  )
  outPutDataFrame=rbind(outPutDataFrame,tempDataFrame)
  
  # Scale:
  MSparameter=matrix(0,2,1)
  
  for (a in 1:2){
    MSparameter[a]=mod.mswm@std[regimenId[a]]
  }
  
  tempDataFrame=data.frame(
    Date=matrix(as.character(tail(Data$Date,1),4,1)),
    Value=MSparameter,
    Ticker=c("stdDevReg1","stdDevReg2")
  )
  outPutDataFrame=rbind(outPutDataFrame,tempDataFrame)
  
  # t smooth prob:
  MSparameter=matrix(0,2,1)
  
  tailSmoProb=tail(mod.mswm@Fit@smoProb,1)
  for (a in 1:2){
    MSparameter[a]=tailSmoProb[regimenId[a]]
  }
  tailSmoProb2=MSparameter
  
  tempDataFrame=data.frame(
    Date=matrix(as.character(tail(Data$Date,1),4,1)),
    Value=MSparameter,
    Ticker=c("smothProbTReg1","smothProbTReg2")
  )
  outPutDataFrame=rbind(outPutDataFrame,tempDataFrame)
  
  # t+1 smooth prob:
  MSparameter=matrix(0,2,1)
  transMat=mod.mswm@transMat
  transMat2=transMat
  for (a in 1:2){
    transMat2[,a]=transMat[,regimenId[a]]
  }
  
  tailSmoProbt1=transMat%*%tailSmoProb2
  tailSmoProbt2=transMat^2%*%tailSmoProb2
  tailSmoProbt3=transMat^3%*%tailSmoProb2
  tailSmoProbt4=transMat^4%*%tailSmoProb2
  tailSmoProbt5=transMat^5%*%tailSmoProb2
  
  tempDataFrame=data.frame(
    Date=matrix(as.character(tail(Data$Date,1),12,1)),
    Value=c(tailSmoProb2,tailSmoProbt1,tailSmoProbt2,tailSmoProbt3,
            tailSmoProbt4,tailSmoProbt5,
            transMat[1,1],transMat[2,1],transMat[1,2],transMat[2,2],
            AIC(mod.mswm),mod.mswm@Fit@logLikel,mod.mswm@k,AIC(mod),logLik(mod),1),
    Ticker=c("smothProbTReg1","smothProbTReg2",
             "smothProbT+1Reg1","smothProbT+1Reg2",
             "smothProbT+2Reg1","smothProbT+2Reg2",
             "smothProbT+3Reg1","smothProbT+3TReg2",
             "smothProbT+4Reg1","smothProbT+4Reg2",
             "smothProbT+5Reg1","smothProbT+5Reg2",
             "transProbS1S1","transProbS1S2","transProbS2S1","transProbS2S2",
             "AIC-MS","LLF-MS","k-MS","AIC-1Reg","LLF-1Reg","k-1Reg")
  )
  
  outPutDataFrame=rbind(outPutDataFrame,tempDataFrame)
  
  # Simple return scenario labeling:
  if (tail(Data$Return,1)>0){
    
    if(tailSmoProb2[2]<=0.5){
      returnScenario2Reg="upside normal" 
      returnScenario2RegN=1    
    } else {
      returnScenario2Reg="upside distress" 
      returnScenario2RegN=2   
    }
    
  } else {
    
    if(tailSmoProb2[2]<=0.5){
      returnScenario2Reg="downside normal" 
      returnScenario2RegN=3  
    } else {
      returnScenario2Reg="downside distress" 
      returnScenario2RegN=4  
    }
    
  }
  
  # single regime return scenario labeling:
  if (tail(Data$Return,1)>0){
    returnScenario1="upside" 
    returnScenario1N=1
  } else {
    returnScenario1="downside"
    returnScenario1N=2}
  
  tempDataFrame=data.frame(Date=c(as.character(tail(Data$Date,1)),
                                  as.character(tail(Data$Date,1))),
                           Value=c(returnScenario1,returnScenario1N),
                           Ticker=c("1RegimeReturnScenario","RegimeReturnScenarioCode"))
  outPutDataFrame=rbind(outPutDataFrame,tempDataFrame)
  
  EstartTime=Sys.time()
  tElaptime=round((EstartTime-TstartTime)/60,4)
  
  # Two regime return scenario labeling:
  tempDataFrame=data.frame(Date=matrix(as.character(tail(Data$Date,1),4,1)),
                           Value=c(returnScenario2Reg,returnScenario2RegN,
                                   as.numeric(estElaptTime),as.numeric(tElaptime)),
                           Ticker=c("twoRegimesReturnScenario","twoRegimesReturnScenarioCode",
                                    "Estimation time","Full date simulation time"))
  outPutDataFrame=rbind(outPutDataFrame,tempDataFrame)  
  
  # Smoothed transition probabilities table:
  smoothTransProb=Data
  
  for (a in 1:2){
    execText=paste0("smoothTransProb$smoothProbReg",a,"=",
           "mod.mswm@Fit@smoProb[1:nrow(Data),regimenId[",a,"]]*100")
    eval(parse(text=execText))
  }
  
  #smoothTransProb$smothProbS1=mod.mswm@Fit@smoProb
    
  outOjb=list(testTable=outPutDataFrame,
              msObject=mod.mswm,
              regmodel=mod,
              outData=smoothTransProb)
  
  cat("\f")
  print(paste0("End of calculations. Ellapsed time: ",as.numeric(tElaptime)," minutes..."))   
  
# Return:
  return(outOjb)
}

