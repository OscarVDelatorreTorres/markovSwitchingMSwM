


mod=lm("Return~1",data=Data)
mod.mswm=msmFit(mod,k=2,p=1,sw=c(T,T,T),control=list(parallel=F))

regimenId=which(mod.mswm@std==max(mod.mswm@std))