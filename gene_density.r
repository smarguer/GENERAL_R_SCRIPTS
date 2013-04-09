#gture=read.delim("splicing_dwn_meiosis_nature.txt",stringsAsFactors=F)
#down_nature=down_nature[,1]
#up_nature=read.delim("splicing_up_meiosis_nature.txt",stringsAsFactors=F)
#up_nature=up_nature[,1]
genes.density<-function(gff, li, win=20, plot=T, write.file=F, file.name="genes.density.output", cor.method="bonferroni",plot.raw=FALSE)

{
	out=data.frame()
	gff=gff[which(gff$"feature"=="gene"),]
        ALL=nrow(gff)
        k=0
	for (i in 1:3)
	{
		print (i)
		wgff=gff[which(gff$"chr"==i),]
		wgff=wgff[order(wgff$"start"),]
		##define matrix##
		for (j in 1:(nrow(wgff)-win))
		{
			k=k+1
			position=wgff[(j+(win/2)),"start"]
			hit = length(which(li %in% wgff$"Name"[j:(j+win-1)]))
			#print(hit)
			test=matrix(0,2,2)
                        test[1,1]=hit
			test[2,1]=win-hit
			#if (test[2,1] < 0)
                        #{
                        # test[2,1]=0
                        #}
                        test[1,2]=length(li)-hit
			test[2,2]=ALL-test[1,1]-test[2,1]-test[1,2]
			fisher = fisher.test(test)
			out[k,1] = i
			out[k,2] = position
			out[k,3] = fisher$p.value
		}
	}

        out[,4]=-log10(out[,3])
        out[,3] = p.adjust(p=out[,3], method=cor.method)
        out[,3] = -log10(out[,3])
        names(out)=c("chr","window.midpoint","cor.p.val","raw.p.val")
        if (plot.raw==FALSE){toplot=3}
        if (plot.raw==TRUE){toplot=4}

        if (plot==T)
        {
	 par(mfrow=c(3,1))
         plot((out[which(out[,1]==1),2]/1000),out[which(out[,1]==1),toplot], ylim=c(0,15), type="l", xlab="position on chromosome [kb]",ylab="-log10 of P value", main="Chromosome 1")
         abline(h=c(1.3), col="orange")
         abline(h=c(3), col="red")
         plot((out[which(out[,1]==2),2]/1000),out[which(out[,1]==2),toplot], ylim=c(0,15), type="l", xlab="position on chromosome [kb]",ylab="-log10 of P value", main="Chromosome 2")
         abline(h=c(1.3), col="orange")
         abline(h=c(3), col="red")
         plot((out[which(out[,1]==3),2]/1000),out[which(out[,1]==3),toplot], ylim=c(0,15), type="l", xlab="position on chromosome [kb]",ylab="-log10 of P value", main="Chromosome 3")
         abline(h=c(1.3), col="orange")
         abline(h=c(3), col="red")
	}
        
        if (write.file==T)
        {
         pdf(file=paste(file.name,"pdf",sep="."))
         
         par(mfrow=c(3,1))
         plot((out[which(out[,1]==1),2]/1000),out[which(out[,1]==1),toplot], ylim=c(0,15), type="l", xlab="position on chromosome [kb]",ylab="-log10 of P value", main="Chromosome 1")
         abline(h=c(1.3), col="orange")
         abline(h=c(3), col="red")
         plot((out[which(out[,1]==2),2]/1000),out[which(out[,1]==2),toplot], ylim=c(0,15), type="l", xlab="position on chromosome [kb]",ylab="-log10 of P value", main="Chromosome 2")
         abline(h=c(1.3), col="orange")
         abline(h=c(3), col="red")
         plot((out[which(out[,1]==3),2]/1000),out[which(out[,1]==3),toplot], ylim=c(0,15), type="l", xlab="position on chromosome [kb]",ylab="-log10 of P value", main="Chromosome 3")
         abline(h=c(1.3), col="orange")
         abline(h=c(3), col="red")
         
         dev.off()

         write.table(out,file=paste(file.name,"txt",sep="."),sep="\t",quote=F,row.names=F)         
        }
        return(out)
}
