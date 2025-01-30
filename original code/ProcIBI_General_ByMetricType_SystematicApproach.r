
# USE THIS IN PLACE OF LINES 137-139 
# for(k in 1:nsamp) {
#   met.samp[k,]<-tapply(metList$METRIC,metList$TYPE,function(x) sample(x,size=1))
# }
# Function that is based on number of metric types in metList
metTypes <- unique(metList$TYPE)

numTypes <- length(metTypes)
expandText <- ''
for(i in 1:length(metTypes)){
  v <- paste0("metList[metList$TYPE==metTypes[", i, "], 'METRIC']")
  if(i==1){
    expandText <- paste0("expand.grid(", v)
  }else{
    expandText <- paste0(expandText, ",", v)
  }
  if(i==length(metTypes)){
    expandText <- paste0(expandText, ")")
  }
}  

met.samp <- eval(parse(text = expandText))
