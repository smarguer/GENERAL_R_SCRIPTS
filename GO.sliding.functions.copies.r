
sliding.GO=function(dat,li,col,win=300,step=1,all=nrow(dat),adjust="fdr")
{
 dat=dat[order(dat[,col]),]
 dat=dat[which(is.na(dat[,col])==FALSE),]
 out=data.frame()
 for(i in 1:(nrow(dat)-win))
 {
  out[i,1]=dat[i,col]
  out[i,2]=dat[i+win,col]
  out[i,3]=dat[i+(round((win/2),0)),col]
  out[i,4]=VENN(dat[i:(i+win),1],li,verbose=FALSE,ALL=all)$p.value 
  out[i,5]=length(which((li %in% dat[i:(i+win),1] == TRUE)))
  out[i,6]=length(which((li %in% dat[i:(i+win),1] == TRUE)))/length(li)*100
  out[i,7]=length(which((li %in% dat[i:(i+win),1] == TRUE)))/win*100
  out[i,8]=VENN(dat[i:(i+win),1],li,verbose=FALSE,ALL=all,alternative="l")$p.value 
  out[i,9]=length(which((dat[i:(i+win),1] %in% li == TRUE)))/length(li)*100
  out[i,10]=length(which((dat[i:(i+win),1] %in% li == TRUE)))/win*100
 }
 out[,4]=p.adjust(out[,4],method=adjust)
 out[,4]=-log10(out[,4])
 out[,8]=p.adjust(out[,8],method=adjust)
 out[,8]=-log10(out[,8])
  
 return(out)
}
#################################################################################
plot.sliding.GO=function(dat,toPlotX=3,toPlotY=7,xlim=c(-5,7),ylim1=c(0,100),ylim2=c(0,100))
{
plot(log2(dat[,toPlotX]),dat[,toPlotY],xlim=xlim,ylim=ylim1, type="l",lwd=3,xlab=expression(bold(paste(Log[2]," mRNA copies/cell"))),ylab=expression(bold("Percent in window")),bty="n",cex.axis=1.35,cex.lab=1.4)
abline(v=0,col="red",lwd=3)
abline(v=c(-1,1),col="red",lwd=2,lty=2)
par(new=T)
plot(log2(dat[,toPlotX]),dat[,4],axes=F,xlab='',ylab='',xlim=xlim,ylim=ylim2,type="l",col=colN,lwd=2,bty="n", cex.axis=1.35,cex.lab=1.4,lty=1)
axis(4,cex.axis=1.35,col=colN,col.axis=colN, line=-2)
mtext(expression(bold(paste(-Log[10]," ",P[Fisher]))),4,line=1,cex=1.2,col=colN)
#lines(log2(dat[,toPlotX]),dat[,8],type="l",col="blue",lwd=2,lty=1)
abline(h=-log10(0.05),lty=2,lwd=1.5)
par(new=T)
plot(log2(dat[,toPlotX]),dat[,toPlotY],xlim=xlim,ylim=ylim1, type="l",lwd=3,bty="n",axes=F,xlab='',ylab='')
abline(v=0,col="red",lwd=3)
abline(v=c(-1,1),col="red",lwd=2,lty=2)

#legend(x=5,y=100,bty="n",legend=length(row.names(GO)[which(GO[,column]==1)]),cex=2)
}



#################################################################################
sliding.GO.screen=function(dat,go,col,win=500,step=100,all=5110,format=T,adjust="fdr")
{
 dat=dat[order(dat[,col]),]
 out=data.frame()
 l=0
 for(i in seq(1,(nrow(dat)-win),step))
 {
  l=l+1
  k=0
  print(i)
  out[l,1]=dat[i,col]
  out[l,2]=dat[i+win,col]
  out[l,3]=dat[i+(round((win/2),0)),col]
  
  for(j in 1:ncol(go))
  {
   li=row.names(go)[which(go[,j]==1)]
   out[l,j+3+k]=VENN(dat[i:(i+win),1],li,verbose=FALSE,ALL=all)$p.value 
   out[l,j+4+k]=length(which((li %in% dat[i:(i+win),1] == TRUE)))
   out[l,j+5+k]=length(which((li %in% dat[i:(i+win),1] == TRUE)))/length(li)*100
   colnames(out)[j+3+k]=paste("P_",colnames(go)[j],sep='')
   colnames(out)[j+4+k]=paste("N_",colnames(go)[j],sep='')
   colnames(out)[j+5+k]=paste("C_",colnames(go)[j],sep='')

   k=k+2
  }
  out[l,grep('^P_',colnames(out))]=p.adjust(out[l,grep('^P_',colnames(out))],method=adjust)
 }
 if(format==TRUE)
 {
  out=out[which(is.na(out[,2])==F),]
  cor=(0.00001*seq(1,nrow(out),1))
  row.names(out)=out[,2]+cor
  out=out[,4:ncol(out)]
  out=t(out)
 }
 #colnames(out)=c("win.start","win.end","win.mid",colnames(go))
 #out[,4]=p.adjust(out[,4],method="holm")
 #out[,4]=-liog10(out[,4])
  
 return(out)
}

##################################################################################

onoD.GO=function(li,dat=BIG,col="MM_prot",go,spread=10,genes=F,tree1=key,tree2=KEY,ext=TRUE,export=FALSE,annot=NULL,plot=TRUE,Quies=FALSE)
{
 if(is.null(annot)==FALSE)
 {
  names=BIG[grep(pattern=annot,BIG[,"Annot"]),1]
 }
 
 else if(length(grep(li,colnames(go)))==0)
 {
  stop("not in GO list")
 }
 
 else if(grepl('GO',li)==TRUE)
 {
  print("Genes assigned to the following terms are shown")
  
  set=c(li,tree1[grep(li,tree1[,"is_a"]),1])
  if(length(set) > 1)
  {
   test.loop=FALSE
   prev=0
   while(test.loop==FALSE)
   {
    cat("#.")
    hold=set
    for(i in 1:length(hold))
    {
     
     hold=unique(c(hold,tree1[grep(hold[i],tree1[,"is_a"]),1]))
     #print(hold)
    }
    set=hold
    if(all(set %in% prev))
    {
     test.loop=TRUE
    }
    else
    {
     prev=set
    }
   }
  }
  cat("\n")
  set1=tree2[which(tree2[,2] %in% set),5]
  
  print(tree2[which(tree2[,2] %in% set),1:4])
  if(length(set1) > 1)
  {
   names=unique(row.names(go)[which(rowSums(go[,set1],na.rm=T)!=0)])
  }
  else
  {
   names=unique(row.names(go)[which(go[,set1]!=0)])
  }
 }
 
 else if(length(ncol(go[,grep(li,colnames(go))]))!=0)
 {
  names=unique(row.names(go)[which(rowSums(go[,grep(li,colnames(go))],na.rm=T)!=0)])
  print(colnames(go)[grep(li,colnames(go))])
 }
 
 else
 {
  names=unique(row.names(go)[which(go[,grep(li,colnames(go))]!=0)])
  print(colnames(go)[grep(li,colnames(go))])
 }
 
 if(plot==TRUE)
 {
  if(Quies==TRUE)
  {
   oneDplot(list(dat[which(dat[,1] %in% names),col],dat[which(dat[,1] %in% names),"MN_prot"]),spread=spread,col="black",ylim=c(5,21),main='')
  }
  else
  {
   oneDplot(list(dat[which(dat[,1] %in% names),col]),spread=spread,col="black",ylim=c(5,21),main='')
  }
  abline(h=c(log2(2300)),col="orange")
  abline(h=c(log2(5000),log2(10000)),col="red")
  abline(h=c(log2(100000),log2(200000)),col="dark green")
  abline(h=log2(40000),col="green")
  abline(h=log2(60000000),col="blue")
  abline(h=log2(100),col="blue")
 }
  
 VIS=c(34,7,8,14,20,21,22,23,33) 
 VIS1=c(1,34,7,8,14,20,21,22,23,33) 
 
 out=BIG[which(BIG[,1] %in% names),VIS]
 out=out[order(out[,col],decreasing=T),]
 out1=BIG[which(BIG[,1] %in% names),VIS1]
 out1=out1[order(out1[,col],decreasing=T),]
 
 if(genes==TRUE)
 {
  print(out)
 }
 
 if(export==TRUE)
 {
  return(out1)
 }
 
 if(ext==F)
 {
  print(BIG[which(BIG[,col]==max(BIG[which(BIG[,1] %in% names),col],na.rm=T)),VIS])
  print(BIG[which(BIG[,col]==min(BIG[which(BIG[,1] %in% names),col],na.rm=T)),VIS])
 }
}

##################################################################################
plot.sliding.average=function(x,y,win=200,step=10,xli=c(2,22),yli=c(2,22),xla='',yla='',med=TRUE,log="xy",main='',cex.leg=2,slideX=c(1:length(x)),add=F)
{
 winh=round((win/2),0)
 test=cbind(x,y)
 test=test[order(test[,1]),]
 slide=vector()
 for(i in 1:nrow(test))
 {
#  print(i)
  if((i > (winh))&((i+(winh)) < (length(test[,1])-winh)))
  {
   slide[i]=median(test[(i-winh):(i+winh),2],na.rm=TRUE)
  }
  else
  {
   slide[i]=NA
  }
 }
 
 if(length(grep('x',log))> 0)
 {
  x[which(x<=0)]=NA
  x=log2(x)
  test[,1]=log2(test[,1])
 }
 if(length(grep('y',log))>0)
 {
  y[which(y<=0)]=NA
  y=log2(y)
  slide=log2(slide)
 }
 
 if(add==T)
 {
 lines(x,y,col="dark grey",pch=20,type="p")
 
 if(med==TRUE)
 {
  lines(test[slideX,1],slide[slideX],lwd=4)
 }
 c=round(cor(x,y,use="pairwise.complete.obs",method="pearson"),2)
 }
 
 else
 { 
 plot(x,y,col="dark grey",pch=20,xlab=xla,ylab=yla,xlim=xli,yli=yli,cex.axis=1.6,cex.lab=2)
 if(med==TRUE)
 {
  lines(test[slideX,1],slide[slideX],lwd=4)
 }
 c=round(cor(x,y,use="pairwise.complete.obs",method="pearson"),2)
 }
 #legend(x=xli[2]/4,y=yli[1]+((yli[2]-yli[1])/5),legend=bquote(cor[pearson] == .(c)),bty="n",cex=cex.leg)
}
##################################
GO.median.table=function(dat,col,go)
{
 out=data.frame()
 for(i in 1:ncol(go))
 {
  li=row.names(go)[which(go[,i]==1)]
  out[i,1]=colnames(go)[i]
  out[i,2]=median(dat[which(dat[,1] %in% li),col],na.rm=TRUE)
  out[i,3]=sd(dat[which(dat[,1] %in% li),col],na.rm=TRUE)
 }
 return(out)
}

##################################################
GOanalysis=function(li,go,all,adjust="fdr",sort=TRUE)
{


 out=data.frame()
 for(i in 1:ncol(go))
 {
  golist=row.names(go)[which(go[,i]==1)]
  out[i,1]=colnames(go)[i]
  out[i,2]=VENN(l2=li,l1=golist,ALL=all,verbose=FALSE)$p.value
  out[i,3]=length(li)
  out[i,4]=length(golist)
  out[i,5]=length(which(li %in% golist))

  #print(i)
 }
 out[,2]=p.adjust(out[,2],method=adjust)
 if(sort==TRUE)
 {
 out=out[order(out[,2]),]
 }
 return(out)
}
#######################################################
sliding.GO.table.copies=function(dat=BIG,col="MM_prot",go=GOtable,win=1000,step=500,adjust="none")
{
 for (i in seq(0,20000,step))
 {
  print(i)
  test=GOanalysis(li=dat[which((dat[,col] > i)&(dat[,col] <= i+win)&(is.na(dat[,col])==FALSE)),1],go=go,all=5110,adjust=adjust,sort=FALSE)
  if(i==0)
  {
   out=data.frame(test[,2])
   rownames(out)=test[,1]
  }
  else
  {
   out=cbind(out,test[,2])
  }
 }
 colnames(out)=seq(0,20000,step)
 return(out)
}

##########################################################
missing=function()
{
MISSING=data.frame()
l=0
for(i in seq(0,5,0.2))
{
print(i)
j=i+0.2
l=l+1
MISSING[l,1]=i
MISSING[l,2]=length(which(is.na(BIG[CODING[which((BIG[CODING,"MM_ex"] >= i)&(BIG[CODING,"MM_ex"] < j))],"MM_prot"])))
MISSING[l,3]=length(which((BIG[CODING,"MM_ex"] >= i)&(BIG[CODING,"MM_ex"] < j)))
MISSING[l,4]=MISSING[l,2]/MISSING[l,3]*100
}
for(i in seq(5,50,5))
{
print(i)
j=i+5
l=l+1
MISSING[l,1]=i
MISSING[l,2]=length(which(is.na(BIG[CODING[which((BIG[CODING,"MM_ex"] >= i)&(BIG[CODING,"MM_ex"] < j))],"MM_prot"])))
MISSING[l,3]=length(which((BIG[CODING,"MM_ex"] >= i)&(BIG[CODING,"MM_ex"] < j)))
MISSING[l,4]=MISSING[l,2]/MISSING[l,3]*100
}
for(i in seq(50,1000,10))
{
print(i)
j=i+10
l=l+1
MISSING[l,1]=i
MISSING[l,2]=length(which(is.na(BIG[CODING[which((BIG[CODING,"MM_ex"] >= i)&(BIG[CODING,"MM_ex"] < j))],"MM_prot"])))
MISSING[l,3]=length(which((BIG[CODING,"MM_ex"] >= i)&(BIG[CODING,"MM_ex"] < j)))
MISSING[l,4]=MISSING[l,2]/MISSING[l,3]*100
}


return(MISSING)
}
############################################################
missing.2=function()
{
MISSING=data.frame()
l=0
for(i in seq(0,5,0.2))
{
print(i)
j=i+0.2
l=l+1
MISSING[l,1]=i
MISSING[l,2]=length(which(is.na(BIG[CODING[which((BIG[CODING,"MM_ex"] >= i)&(BIG[CODING,"MM_ex"] < j))],"MM_prot"])))
MISSING[l,3]=length(which((BIG[CODING,"MM_ex"] >= i)&(BIG[CODING,"MM_ex"] < j)))
MISSING[l,4]=MISSING[l,2]/MISSING[l,3]*100
}
for(i in seq(5,50,5))
{
print(i)
j=i+5
l=l+1
MISSING[l,1]=i
MISSING[l,2]=length(which(is.na(BIG[CODING[which((BIG[CODING,"MM_ex"] >= i)&(BIG[CODING,"MM_ex"] < j))],"MM_prot"])))
MISSING[l,3]=length(which((BIG[CODING,"MM_ex"] >= i)&(BIG[CODING,"MM_ex"] < j)))
MISSING[l,4]=MISSING[l,2]/MISSING[l,3]*100
}
for(i in seq(50,1000,10))
{
print(i)
j=i+10
l=l+1
MISSING[l,1]=i
MISSING[l,2]=length(which(is.na(BIG[CODING[which((BIG[CODING,"MM_ex"] >= i)&(BIG[CODING,"MM_ex"] < j))],"MM_prot"])))
MISSING[l,3]=length(which((BIG[CODING,"MM_ex"] >= i)&(BIG[CODING,"MM_ex"] < j)))
MISSING[l,4]=MISSING[l,2]/MISSING[l,3]*100
}


return(MISSING)
}
############################################################
missing.3=function(go=GOfinal)
{
 MISSING=data.frame()
 for(i in 1:ncol(go))
 {
  print(i)

  MISSING[i,1]=length(which(is.na(BIG[CODING[which(BIG[CODING,1] %in% row.names(go)[which(go[,i] ==1)])],"MM_prot"])==TRUE))
  MISSING[i,2]=length(row.names(go)[which(go[,i] ==1)])
  MISSING[i,3]=MISSING[i,1]/MISSING[i,2]*100
  MISSING[i,4]=100-MISSING[i,3]
  row.names(MISSING)[i]=colnames(go)[i]
 }
 return(MISSING)
}

#############################################################
missing.4=function(go=GOfinal)
{
 MISSING=data.frame()
 for(i in 1:ncol(go))
 {
  print(i)
  tot.low=length(BIG[CODING[which((BIG[CODING,"MM_ex"] >= 0)&(BIG[CODING,"MM_ex"] < 0.5)&(BIG[CODING,1] %in% row.names(go)[which(go[,i] ==1)]))],"MM_prot"])
  MISSING[i,1]=length(which(is.na(BIG[CODING[which((BIG[CODING,"MM_ex"] >= 0)&(BIG[CODING,"MM_ex"] < 0.5)&(BIG[CODING,1] %in% row.names(go)[which(go[,i] ==1)]))],"MM_prot"])==TRUE))
  MISSING[i,2]=tot.low-MISSING[i,1]
  
  tot.mid=length(BIG[CODING[which((BIG[CODING,"MM_ex"] >= 0.5)&(BIG[CODING,"MM_ex"] < 2)&(BIG[CODING,1] %in% row.names(go)[which(go[,i] ==1)]))],"MM_prot"])
  MISSING[i,3]=length(which(is.na(BIG[CODING[which((BIG[CODING,"MM_ex"] >= 0.5)&(BIG[CODING,"MM_ex"] < 2)&(BIG[CODING,1] %in% row.names(go)[which(go[,i] ==1)]))],"MM_prot"])==TRUE))
  MISSING[i,4]=tot.mid-MISSING[i,3]

  tot.high=length(BIG[CODING[which((BIG[CODING,"MM_ex"] >= 2)&(BIG[CODING,"MM_ex"] < 1000)&(BIG[CODING,1] %in% row.names(go)[which(go[,i] ==1)]))],"MM_prot"])

  MISSING[i,5]=length(which(is.na(BIG[CODING[which((BIG[CODING,"MM_ex"] >= 2)&(BIG[CODING,"MM_ex"] < 1000)&(BIG[CODING,1] %in% row.names(go)[which(go[,i] ==1)]))],"MM_prot"])==TRUE))
  MISSING[i,6]=tot.high-MISSING[i,5]

  row.names(MISSING)[i]=colnames(go)[i]
  
 }
 return(MISSING)
}


###########################################################
cumul.plot=function(dat=BIG,col="MM")
{
 col1=paste(col,"_prot",sep='')
 col2=paste(col,"_ex",sep='')
 dat1=cbind(dat[,col1],dat[,col2])
 dat1=dat1[order(dat1[,1],decreasing=T),]
 dat1=dat1[which(is.na(dat1[,1])==F),]
 #dat2=dat2[which(is.na(dat1)==F)]
 out=data.frame(stringsAsFactors=T)
 for(i in 1:nrow(dat1))
 {
  out[i,1]=i/nrow(dat1)*100
  out[i,2]=sum(dat1[1:i,1],na.rm=T)/sum(dat1[,1],na.rm=T)*100
  out[i,3]=sum(dat1[1:i,2],na.rm=T)/sum(dat[,col2],na.rm=T)*100
 }
  #par(mfrow=c(1,2))
  print(min(which(out[,2]>=50)))
  plot(out[,1],out[,3],cex.axis=2,cex.lab=2,xlab=expression(bold("Protein abundance [% rank]")), ylab=expression(bold("Cumulative copies [%]")),type="l",lwd=4,col="black", ylim=c(0,100))
  lines(out[,1],out[,2],type="l",lwd=7,col="red")
  abline(v=c(20),lwd=3,lty=2,col="blue")
  abline(h=c(80),lwd=3,lty=2,col="green")
  legend(x="center",legend=c(expression(bold("Proteins")),expression(bold("mRNAs"))),text.col=c("red","black"),cex=2,bty="n")

  return(out)
}


###########################################################
cumul.plot1=function(dat=BIG,col="MM")
{
 #col1=paste(col,"_prot",sep='')
 col2=paste(col,"_ex",sep='')
 #dat1=cbind(dat[,col1],dat[,col2])
 dat1=dat[,col2]
 dat1=dat1[order(dat1,decreasing=T)]
 dat1=dat1[which(is.na(dat1)==F)]
 #dat2=dat2[which(is.na(dat1)==F)]
 out=data.frame(stringsAsFactors=T)
 for(i in 1:length(dat1))
 {
  out[i,1]=i/length(dat1)*100
  out[i,2]=sum(dat1[1:i],na.rm=T)/sum(dat1,na.rm=T)*100
  #out[i,3]=sum(dat1[1:i,2],na.rm=T)/sum(dat[,col2],na.rm=T)*100
 }
  #par(mfrow=c(1,2))
  print(min(which(out[,2]>=50)))
  plot(out[,1],out[,2],cex.axis=2,cex.lab=2,xlab=expression(bold("lncRNA abundance [% rank]")), ylab=expression(bold("Cumulative lncRNA copies [%]")),type="l",lwd=4,col="black", ylim=c(0,100))
  #lines(out[,1],out[,2],type="l",lwd=7,col="red")
  abline(v=c(10,50,90),lwd=1,lty=1,col="blue")
  abline(h=c(10,50,90),lwd=1,lty=1,col="green")
  abline(v=c(14.5),lwd=3,lty=2,col="red")
  #plot(out[,1],log2(dat1[,1]),cex.axis=1.6,cex.lab=1.6,xlab="percent of ranked detected proteins", ylab="Log2 proteins copies",type="l",lwd=7,col="red")
  #plot(out[1:200,2])
  return(out)
}

###########################################################
cumul.plot2=function(dat=BIG[CODING,])
{
 col1="MM_prot"
 col2="MN_prot"
 dat1=cbind(dat[,col1],dat[,col2])
 #dat1=dat[,col2]
 dat1=dat1[order(dat1[,1],decreasing=T),]
 dat1=dat1[which(is.na(dat1[,1])==F),]
 dat1[which(is.na(dat1[,2])==T),2]=0
 #dat2=dat2[which(is.na(dat1)==F)]
 out=data.frame(stringsAsFactors=T)
 for(i in 1:nrow(dat1))
 {
  out[i,1]=i/nrow(dat1)*100
  #if(dat1[i,1]-dat1[i,2] < 0)
  #{
   out[i,2]=sum(dat1[1:i,1]-dat1[1:i,2],na.rm=T)/(sum(dat1[,1])-sum(dat1[,2]))*100
  #}
  #else if (i==1)
  #{
  # out[i,2]=0
  #}
  #else
  #{
  # out[i,2]=out[i-1,2]
  #}
  #out[i,3]=sum(dat1[1:i,2],na.rm=T)/sum(dat[,col2],na.rm=T)*100
 }
  #par(mfrow=c(1,2))
  print(min(which(out[,2]>=50)))
  plot(out[,1],out[,2],cex.axis=1.6,cex.lab=1.6,xlab="percent of ranked proteins during prolifeartion", ylab="number of protein lost",type="l",lwd=4,col="black")
  #lines(out[,1],out[,2],type="l",lwd=7,col="red")
  abline(v=c(10,50,90),lwd=1,lty=1,col="blue")
  abline(h=c(10,50,90),lwd=1,lty=1,col="green")
  #abline(v=c(14.5),lwd=3,lty=2,col="red")
  #plot(out[,1],log2(dat1[,1]),cex.axis=1.6,cex.lab=1.6,xlab="percent of ranked detected proteins", ylab="Log2 proteins copies",type="l",lwd=7,col="red")
  #plot(out[1:200,2])
  return(out)
}



############################################################

Dplot_copies=function(f1, p=0.5, ncolours=21, li=T, leg=TRUE, xl=colnames(f1)[1],yl=colnames(f1)[2], m="plot", scale=c(-10,15),cex.axis=1.6)
{
 par(bg="white", fg="black")
 plot(log2(f1), xlim=scale, ylim=scale, type="n", col.lab="black", col.axis="black", cex.lab=2, xlab=xl,ylab=yl, main=m,cex.axis=cex.axis,pch=20)
 Image(log2(f1),colramp=topo.colors,  pixs=p)

 if(li==T)
 {
  abline(0,1,lwd=1)
  abline(1,1,lwd=1,lty=2)
  abline(-1,1,lwd=1,lty=2)
  abline(v=0,h=0,col="red",lwd=1.5)
  abline(v=c(-1,1),h=c(-1,1),col="red",lwd=1,lty=2)
  #abline(v=c(log2(1),log2(10),log2(100)), col="orange")
  #abline(h=c(log2(1),log2(10),log2(100)), col="orange")
 }

 if(leg==TRUE)
 {
  legend(x=-10,y=14,legend=paste("cor=",signif(cor(log2(f1[,1]),log2(f1[,2]), use="pairwise"),2),sep=''), bty="n", cex=1.5)
  #legend(x=-10,y=14,legend=paste("cor=",signif(cor((f1[,1]),(f1[,2]), use="pairwise"),2),sep=''), bty="n", cex=1.5)
 } 
}

Dplot_copies_2=function(f1, p=0.5, ncolours=21, li=T, leg=TRUE, xl=colnames(f1)[1],yl=colnames(f1)[2], m="plot", xli=c(-5,15),yli=c(5,25),cex.axis=1.6,cex.lab=2)
{
 par(bg="white", fg="black")
 plot(log2(f1), xlim=xli, ylim=yli, type="n", col.lab="black", col.axis="black", cex.lab=cex.lab,cex.axis=cex.axis, xlab=xl,ylab=yl, main=m)
 Image(log2(f1),colramp=topo.colors,  pixs=p)

 if(li==T)
 {
  abline(0,1,lwd=1)
  abline(1,1,lwd=1,lty=2)
  abline(-1,1,lwd=1,lty=2)
  abline(v=c(log2(1),log2(10),log2(100)), col="orange")
  abline(h=c(log2(1),log2(10),log2(100)), col="orange")
 }

 if(leg==TRUE)
 {
  #legend(x=-10,y=14,legend=paste("cor=",signif(cor(log2(f1[,1]),log2(f1[,2]), use="pairwise"),3),sep=''), bty="n", cex=1.5)
  legend(x=-10,y=14,legend=paste("cor=",signif(cor((f1[,1]),(f1[,2]), use="pairwise"),3),sep=''), bty="n", cex=1.5)
 } 
}



















