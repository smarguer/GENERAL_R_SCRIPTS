




signal.along.chromosome=function(chr, size=1000, shift=100)
{
out=matrix(0,((length(chr)/shift)),2)

for (i in 0:((length(chr)/shift)))
{
print(i)
start=i*shift
if (start==0)
{
start=1
}
end=(i*shift)+size
if ((start+size) > length(chr))
{
next
}

out[i,1]=start
out[i,2]=median(chr[start:end])/1000
}

#print(out[,1],out[,2])
return(out)
}


genes.along.chromosome=function(exp, cond, chr, size=10, shift=1, gff=mygff, what="med")
{

gffg=gff[which(gff[,"feature"]=="gene"),]
gffg=gffg[which(gffg[,"chr"]==chr),]

out=matrix(0,((nrow(gffg)/shift)-size),3)

for (i in 1:((nrow(gffg)/shift)-size))
{
#print(i)

start=i*shift+1
if (start==0)
{
start=1
}
end=(i*shift)+size

if ((start+size) > nrow(gff))
{
print("END")
}
out[i,1]=i
out[i,2]=gffg[(i+(size/2)),"start"]/1000
if (what=="med")
{
out[i,3]=median(exp[which(row.names(exp) %in% gffg[start:end,"Name"]),cond], na.rm=T)
}
else if (what=="sum")
{
out[i,3]=sum(exp[which(row.names(exp) %in% gffg[start:end,"Name"]),cond], na.rm=T)
}




}
if (print==T)
{
#jpeg("figure_pab2_020610_3.jpg", height=35, width=25, unit="cm", res=600)
par(mfrow=c(3,1))
plot(genes.wt.1[,2],log2(genes.pab2.1[,3]/genes.wt.1[,3]), type="l", ylim=c(-5,10), xlab="coordinates [kb]", ylab="mutant to wt ratio", main="Chromosome 1")
abline(h=c(-1,1),col="red")
abline(v=c(3753.687,3789.421), col="grey")
plot(genes.wt.2[,2],log2(genes.pab2.2[,3]/genes.wt.2[,3]), type="l", ylim=c(-5,10), xlab="coordinates [kb]", ylab="mutant to wt ratio", main="Chromosome 2")
abline(h=c(-1,1), col="red")
abline(v=c(1604.786,1638.232), col="grey")
plot(genes.wt.3[,2],log2(genes.pab2.3[,3]/genes.wt.3[,3]), type="l", ylim=c(-5,10), xlab="coordinates [kb]", ylab="mutant to wt ratio", main="Chromosome 3")
abline(h=c(-1,1), col="red")
abline(v=c(1073.708,1137.003), col="grey")
dev.off()
}
#print(out[,1],out[,2])
return(out)
}




