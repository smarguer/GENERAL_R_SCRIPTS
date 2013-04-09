VENN=function(l1,l2,ALL=5880, plot=F, verbose=T)
{
if (plot==T)
{
doVennDiagram(a=l1,b=l2)
}
test=matrix(0,2,2)
hit=length(which(l1 %in% l2))
test[1,1]=hit
test[2,1]=length(l1)-hit
test[1,2]=length(l2)-hit
test[2,2]=ALL-test[1,1]-test[2,1]-test[1,2]
fisher = fisher.test(test, alternative="g")
if(verbose==T)
{
print(paste("(",test[2,1],"(",hit,")",test[1,2],")",sep=""))
}

return(fisher)
}

####################
VENN.clusters=function(clusters,L,ALL=5880,plot=F)
{
out=data.frame()

for (i in 1:length(clusters))
{

cluster.test=clusters[[i]][which(clusters[[i]]!="")]

test=matrix(0,2,2)
hit=length(which(cluster.test %in% L))
test[1,1]=hit
test[2,1]=length(cluster.test)-hit
test[1,2]=length(L)-hit
test[2,2]=ALL-test[1,1]-test[2,1]-test[1,2]
fisher = fisher.test(test, alternative="g")
out[i,1]=paste("cluster",i,sep="")
out[i,2]=length(cluster.test)
out[i,3]=length(L)
out[i,4]=hit
out[i,5]=fisher$p.value
}
out[,5]=p.adjust(p=out[,5], method="bonferroni")
if(plot==T)
{
barplot(-log10(out[,5]),ylim=c(0,10),names.arg=seq(1,length(clusters),1))
abline(h=-log10(0.05),col="blue")
abline(h=-log10(0.001),col="red")
}
return(out)
}
