cat(
"###HELP###","\n",
"'oneDplot()' plots 1D dot plots.","\n",
"Takes up f1, a list of vectors to be plotted.","\n",
"Takes up only linear values.","\n",
"Currently only process values in log space. This means that diplaying linear values won't look that good and that values equal to zero are discarded.","\n",
"Three parameter define how the distribution looks like:","\n",
"breaks: similar to function hist. Splits data distribution in 'breaks' number of bins.","\n",
"spread: horizontal spread of the data points (for the moment the lower the number the larger the spread).","\n",
"pwidth: distance between plots.","\n",
"Tweak breaks,spread,pwidth to get good looking plots.","\n",
"'sep' will draw lines between plots.","\n",
sep='')

#"to get several plots next to each other input a table. To overlay plots set 'add' to TRUE.","\n",
#"Takes up only linear values, to display log values use log argument.","\n",
oneDplot=function(f1, ylim=NULL, log=TRUE, add=FALSE, sep=FALSE, spread=10, breaks=15, col="light blue", pch=20,pwidth=5,ylab='',xlab='',main=NULL, names='',cex=1,smooth=NULL,highlight=NULL, highlight.col="red",cex.lab=1,cex.axis=1,cex.labels=1,las=0)
{

#####a few adjustments to deal with differing inputs whether one or more datasets are to be plotted
 if (is.list(f1)==FALSE)
 {
  stop ("input needs to be a list. Write list(your datasets seprated by comas).")
 } 
 toDo=length(f1)
 
 if (length(spread)!=length(f1))
 {
  spread=rep(spread,(toDo/length(spread)))
 }
 
 if (length(names)!=length(f1))
 {
  names=rep(names,(toDo/length(names)))
 }
 if (length(col)!=length(f1))
 {
  col=rep(col,(toDo/length(col)))
 }
 if (length(pch)!=length(f1))
 {
  pch=rep(pch,(toDo/length(pch)))
 }
 if (length(highlight.col)!=length(f1))
 {
  highlight.col=rep(highlight.col,(toDo/length(highlight.col)))
 }
 if (is.null(smooth)==FALSE)
 {
  cat("smoothing distribution...","\n",sep='')
 }

#######defines axis range
 xlim=c(-pwidth,(0+(toDo*pwidth)))

 if(is.null(ylim)==TRUE)
 {
  #ylim=c((mi-(floor(mi/5))),(max+(ceiling(ma/5))))
  ylim=c(-10,30)
 }

#####computes spread for each of the 'toDo' datasets
  
 for (j in 1:toDo)
 {
  data=f1[[j]]
  data1=rep(1,length(data))
  data2=rep(0,length(data))
  data=cbind(data1,data,data2)
  
  if (is.null(highlight)==FALSE)
  {
   data[highlight[[j]],3]=1 
  }
  
  data=data[which(is.na(data[,2])==FALSE),]
  
  if(is.null(nrow(data))==TRUE)
  {
   if(log==TRUE)
   {
    data[2]=log2(data[2])
   }
   if(add==TRUE)
   { 
    lines((0+(pwidth*(j-1))),data[2], col=col[j],pch=pch[j], xaxt="n",cex=cex, type="p")
   }
   else
   {
    plot((0+(pwidth*(j-1))),data[2], ylim=ylim, xlim=xlim, col=col[j],pch=pch[j], xaxt="n",xlab=xlab, ylab=ylab, main=main,cex=cex)
   }
   next
  }
  
  data=data[order(data[,2]),]
  ma=max(data[,2])
  mi=min(data[,2])
 
  range=ceiling(log2(ma-mi+1))
  ma=ceiling(log2(ma))
  mi=floor(log2(mi))  
  
  #####'inter' = boundaries of argument 'breaks' number of bins
  bin=(ma-mi)/(breaks-1) 
  inter=seq(mi,ma,bin)
#print(inter) 
  
  for (i in inter)
  {
   start=i
   end=i+bin
   
   if (is.null(smooth)==FALSE)
   {
    end=i+(smooth*bin)
   }
   
   if(i != ma)
   {
    slice=which((log2(data[,2])>=start) & (log2(data[,2])<end))
   }
   else if (i == ma)
   {
    slice=which((log2(data[,2])>=start) & (log2(data[,2])<=end))
   }

   nu=length(slice)

   if(nu > 1)
   {
    dat=((seq(1,nu,1)-median(seq(1,nu,1),na.rm=TRUE))/spread[j]) + (pwidth * (j-1))
    data[slice,1]=sample(dat,size=length(dat))
   }

   else if (nu >= 0)
   {
    data[slice,1]=0 + (pwidth * (j-1))
   }
  }
 
#print(data)


#######plots data

  if (log==TRUE)
  {
   data[,2]=log2(data[,2])
  }
 
  if((j==1) & (add==FALSE))
  {
   plot(data[,1],data[,2], ylim=ylim, xlim=xlim, col=col[j],pch=pch[j], xaxt="n",xlab=xlab, ylab=ylab, main=main,cex=cex,cex.lab=cex.lab,cex.axis=cex.axis)
  }
 
  else if ((j!=1)|(add==TRUE))
  {
   lines(data[,1],data[,2], col=col[j],pch=pch[j], type="p",cex=cex)
  }
  
  if (is.null(highlight)==FALSE)
  {
   lines(data[which(data[,3]==1),1],data[which(data[,3]==1),2],col=highlight.col[j],pch=pch[j],type="p",cex=cex)
  }
 }

 axis(side=1,at=seq(0,((toDo-1)*pwidth),pwidth),labels=names,cex.axis=cex.labels,las=las)

######draws lines between plot if 'sep' is TRUE
 if (sep==TRUE)
 {
  seps=seq(-(pwidth/2),((toDo*pwidth)+pwidth/2),pwidth)
  abline(v=seps)
 }
}




