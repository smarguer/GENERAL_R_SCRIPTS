plot.geneOrientation=function(output=out,title="Default", ylim1=c(0,10))
{
plot(row.names(output),-log10(output[,1]),ylim=ylim1,type="p",pch=20,lwd=3, xlab="ranked bin", ylab="-log10 binomial p.value", main=title)
lines(row.names(output),-log10(output[,2]),type="p",col="blue",pch=20,lwd=3)
lines(row.names(output),-log10(output[,3]),type="p",col="green",pch=20,lwd=3)
lines(row.names(output),-log10(output[,4]),type="p",col="pink",pch=20,lwd=3)
#lines(row.names(output),-log10(output[,5]),type="p",col="orange",pch=20,lwd=3)
#lines(row.names(output),-log10(output[,6]),type="p",col="red",pch=20,lwd=3)
#lines(row.names(output),-log10(output[,7]),type="p",col="purple",pch=20,lwd=3)
abline(h=c(1.3,3))
#legend(x="topright",legend=colnames(output),fill=c("black","blue","green","pink","orange","red","purple"))
legend(x="topright",legend=colnames(output)[1:4],fill=c("black","blue","green","pink","orange","red","purple"))
}



sliding.geneOrientation=function(data1=test,data2=RPKM,totest="il_wt1",li=c(1:nrow(data2)),sreal=T,win=100, step=50, alt="greater",adjust=T, plot="gl")
{
data2=data2[li,]
list1=data2[order(data2[,totest],decreasing=T),"name"]
outG=matrix(0,7,((nrow(data1)-win)/step)+1)
outL=matrix(0,7,((nrow(data1)-win)/step)+1)
bins=rep(0,((nrow(data1)-win)/step)+1)
i=0
j=0
while (i <=(nrow(data1)-win))
{
k=i
if(k==0){k=1}
i=i+step
j=j+1
#print(paste(j,i,sep=" "))
bin.test=bin.geneOrientation(data=data1,list=which(data1[,1] %in% list1[k:(k+win-1)]),real=sreal,plot=F,verbose=F)
outG[,j]=as.numeric(bin.test[,"greater"])
outL[,j]=as.numeric(bin.test[,"less"])
bins[j]=k
}

outG=t(outG)
outG=as.data.frame(outG)
outL=t(outL)
outL=as.data.frame(outL)

if(adjust==T)
{
 for (l in 1:7)
 {
  outG[,l]=p.adjust(p=outG[,l], method="bonferroni",n=7*nrow(outG))
  outL[,l]=p.adjust(p=outL[,l], method="bonferroni",n=7*nrow(outL))
 }
}
colnames(outG)=c("Divergent[1,PP]:","Convergent[2,TT]:","tandem AS[3,PT]:","tandem S[4,TP]:","tadem [3+4]:","prom/prom[1+3]:","ter/ter[2+3]:")
row.names(outG)=bins
colnames(outL)=c("Divergent[1,PP]:","Convergent[2,TT]:","tandem AS[3,PT]:","tandem S[4,TP]:","tadem [3+4]:","prom/prom[1+3]:","ter/ter[2+3]:")
row.names(outL)=bins
if ((plot=="gl")|(plot=="g"))
{
plot.geneOrientation(output=outG,title=paste(totest,"greater",sep=" "), ylim1=c(0,10))
}
if ((plot=="gl")|(plot=="l"))
{
plot.geneOrientation(output=outL,title=paste(totest,"less",sep=" "), ylim1=c(0,10))
}
return(outG)
}


bin.geneOrientation=function(data,list,real=F,plot=T,verbose=T)
{


observed=c(length(which(data[list,4]==1)),length(which(data[list,4]==2)),length(which(data[list,4]==3)),length(which(data[list,4]==4)))
all=sum(observed)
con=(observed[2])
div=(observed[1])
mix=observed[3]+observed[4]
ter=observed[2]+observed[3]
prom=observed[1]+observed[3]
expected=rep(sum(observed)/4,4)
prob=c(0.25,0.25,0.25,0.25)
if(real==T)
{
prob=c(0.244,0.244,0.279,0.233)
expected=c(nrow(data[list,])*prob[1],nrow(data[list,])*prob[2],nrow(data[list,])*prob[3],nrow(data[list,])*prob[4])
}
tableout=cbind(observed,expected)
row.names(tableout)=c("PP[1]","TT[2]","PT[3]","TP[4]")
colnames(tableout)=c("obs","exp")

if(verbose==T)
{
print(tableout)
}

bin.con.g=binom.test(x=c(con,(all-con)),p=prob[2],alternative="greater")
bin.div.g=binom.test(x=c(div,(all-div)),p=prob[1],alternative="greater")
bin.mix.g=binom.test(x=c(mix,(all-mix)),p=(prob[3]+prob[4]),alternative="greater")
bin.prom.g=binom.test(x=c(prom,(all-prom)),p=(prob[1]+prob[3]),alternative="greater")
bin.ter.g=binom.test(x=c(ter,(all-ter)),p=(prob[2]+prob[3]),alternative="greater")
bin.PT.g=binom.test(x=c(observed[3],(all-observed[3])),p=prob[3],alternative="greater")
bin.TP.g=binom.test(x=c(observed[4],(all-observed[4])),p=prob[4],alternative="greater")

bin.con.l=binom.test(x=c(con,(all-con)),p=prob[2],alternative="less")
bin.div.l=binom.test(x=c(div,(all-div)),p=prob[1],alternative="less")
bin.mix.l=binom.test(x=c(mix,(all-mix)),p=(prob[3]+prob[4]),alternative="less")
bin.prom.l=binom.test(x=c(prom,(all-prom)),p=(prob[1]+prob[3]),alternative="less")
bin.ter.l=binom.test(x=c(ter,(all-ter)),p=(prob[2]+prob[3]),alternative="less")
bin.PT.l=binom.test(x=c(observed[3],(all-observed[3])),p=prob[3],alternative="less")
bin.TP.l=binom.test(x=c(observed[4],(all-observed[4])),p=prob[4],alternative="less")


out=data.frame()
g=c(bin.div.g$p.value,bin.con.g$p.value,bin.PT.g$p.value,bin.TP.g$p.value,bin.mix.g$p.value,bin.prom.g$p.value,bin.ter.g$p.value)
l=c(bin.div.l$p.value,bin.con.l$p.value,bin.PT.l$p.value,bin.TP.l$p.value,bin.mix.l$p.value,bin.prom.l$p.value,bin.ter.l$p.value)
t=c("Divergent[1,PP]:","Convergent[2,TT]:","tandem AS[3,PT]:","tandem S[4,TP]:","tadem [3+4]:","prom/prom[1+3]:","ter/ter[2+3]:")
out=cbind(g,l)
out=as.data.frame(out)
row.names(out)=t
colnames(out)=c("greater","less")
out=format(out,digits=3)

if (plot==T)
{
par(mfrow=c(1,3))
hist(as.numeric(data[list,4]), breaks=100,ylim=c(0,(max(tableout[,1])+max(tableout[,1])/4)))
abline(h=tableout[1,2],col="red")
barplot(height=-log10(as.numeric(out[,"greater"])),beside=T, main="binomial test, alternative GREATER",names.arg=c("div","con","tan.as","tan.s","tan","prom","ter"),ylim=c(0,15))
abline(h=1.3,col="blue")
abline(h=3,col="red")
barplot(height=-log10(as.numeric(out[,"less"])),beside=T, main="binomial test, alternative LESS",names.arg=c("div","con","tan.as","tan.s","tan","prom","ter"),ylim=c(0,15))
abline(h=1.3,col="blue")
abline(h=3,col="red")
}


return(out)
}




chi.geneOrientation=function(data,list)
{
observed=c(length(which(data[list,4]==1)),length(which(data[list,4]==2)),length(which(data[list,4]==3)),length(which(data[list,4]==4)))
all=sum(observed)
con=(observed[2])
div=(observed[1])
mix=observed[3]+observed[4]
ter=observed[2]+observed[3]
prom=observed[1]+observed[3]
expected=rep(sum(observed)/4,4)

tableout=cbind(c(1,2,3,4),observed,expected)
colnames(tableout)=c("cat","obs","exp")
print(tableout)
hist(as.numeric(data[list,4]), breaks=100,ylim=c(0,(max(tableout[,2])+max(tableout[,2])/4)))
#hist(as.numeric(data[list,4]), breaks=100)
abline(h=tableout[1,3],col="red")
#table=cbind(observed,expected)
chi.tot=chisq.test(x=observed,p=c(0.25,0.25,0.25,0.25))
chi.con=chisq.test(x=c(con,(all-con)),p=c(0.25,0.75))
chi.div=chisq.test(x=c(div,(all-div)),p=c(0.25,0.75))
chi.mix=chisq.test(x=c(mix,(all-mix)),p=c(0.5,0.5))
chi.prom=chisq.test(x=c(prom,(all-prom)),p=c(0.5,0.5))
chi.ter=chisq.test(x=c(ter,(all-ter)),p=c(0.5,0.5))
chi.PT=chisq.test(x=c(observed[3],(all-observed[3])),p=c(0.25,0.75))
chi.TP=chisq.test(x=c(observed[4],(all-observed[4])),p=c(0.25,0.75))

print(paste("all:",chi.tot$p.value,sep=' '))
print("****")
print(paste("Divergent[1,PP]:",chi.div$p.value,sep=' '))
print(paste("Convergent[2,TT]:",chi.con$p.value,sep=' '))
print(paste("tandem AS[3,PT]:",chi.PT$p.value,sep=' '))
print(paste("tandem S[4,TP]:",chi.TP$p.value,sep=' '))
print(paste("tadem:",chi.mix$p.value,sep=' '))
print("****")
print(paste("ter/ter[2+3]:",chi.ter$p.value,sep=' '))
print(paste("prom/prom[1+3]:",chi.prom$p.value,sep=' '))
}












geneOrientation=function(gffa=gff)
{

gffg=gffa[which(gffa$feature=="gene"),]
out=data.frame()

for (i in 1:nrow(gffg))
{

j=i-1
k=i+1
if (j<1){j=i}
if (k>nrow(gffg)){k=i}

out[i,1]=gffg$Name[i]
out[i,2]=NA
out[i,3]=NA
out[i,4]=NA
out[i,5]=gffg$chr[i]

if (gffg$strand[i]=="+")
{
 if (gffg$strand[j]=="+") 
 {
  if (gffg$strand[k]=="+")
  {
   out[i,2]="T"
   out[i,3]="P"
   out[i,4]=4
  }
  if (gffg$strand[k]=="-")
  {
   out[i,2]="T"
   out[i,3]="T"
   out[i,4]=2
  }
 }

 if (gffg$strand[j]=="-") 
 {
  if (gffg$strand[k]=="+")
  {
   out[i,2]="P"
   out[i,3]="P"
   out[i,4]=1
  }
  if (gffg$strand[k]=="-")
  {
   out[i,2]="P"
   out[i,3]="T"
   out[i,4]=3
  }
 }
}
if (gffg$strand[i]=="-")
{
 if (gffg$strand[j]=="+") 
 {
  if (gffg$strand[k]=="+")
  {
   out[i,2]="P"
   out[i,3]="T"
   out[i,4]=3
  }
  if (gffg$strand[k]=="-")
  {
   out[i,2]="T"
   out[i,3]="T"
   out[i,4]=1
  }
 }

 if (gffg$strand[j]=="-") 
 {
  if (gffg$strand[k]=="+")
  {
   out[i,2]="P"
   out[i,3]="P"
   out[i,4]=2
  }
  if (gffg$strand[k]=="-")
  {
   out[i,2]="T"
   out[i,3]="P"
   out[i,4]=4
  }
 }
}
}
return(out)
}
