plot.list.along<-function(gff, li, plot=T, write.file=F, file.name="plot.list.along.output")
{
	if(is.null(ncol(li)))
        {
         stop("I need two columns, first systematic gene names, second scores") 
        }
	if(ncol(li) !=2)
        {
         stop("I only need two columns, first systematic gene names, second scores") 
        }


        out=data.frame()
	gff=gff[which(gff$"feature"=="gene"),]
        k=0
	for (i in 1:3)
	{
		print (i)
		wgff=gff[which(gff$"chr"==i),]
		wgff=wgff[order(wgff$"start"),]
		##define matrix##
                for (j in 1:nrow(wgff))
		{
                        k=k+1	
			position=wgff[j,"start"]
	                if (wgff[j,"Name"] %in% li[,1])
			{
                         score=li[which(li[,1]==wgff[j,"Name"]),2]
                        }
			else
                        {
                         score=0
                        }
			out[k,1] = i
			out[k,2] = position
			out[k,3] = score
		}
	}

        names(out)=c("chr","ORF.start","score")
        if (plot==T)
        {
	 par(mfrow=c(3,1))
         plot((out[which(out[,1]==1),2]/1000),out[which(out[,1]==1),3], ylim=c(-10,10), type="h", xlab="position on chromosome [kb]",ylab="score", main="Chromosome 1")
         abline(h=c(0), col="black")
         abline(h=c(-2,2), col="red")
         abline(v=c(3753.687,3789.421), col="grey")
         plot((out[which(out[,1]==2),2]/1000),out[which(out[,1]==2),3], ylim=c(-10,10), type="h", xlab="position on chromosome [kb]",ylab="score", main="Chromosome 2")
         abline(h=c(0), col="black")
         abline(h=c(-2,2), col="red")
         abline(v=c(1604.786,1638.232), col="grey")
         plot((out[which(out[,1]==3),2]/1000),out[which(out[,1]==3),3], ylim=c(-10,10), type="h", xlab="position on chromosome [kb]",ylab="score", main="Chromosome 3")
         abline(h=c(0), col="black")
         abline(h=c(-2,2), col="red")
         abline(v=c(1073.708,1137.003), col="grey")
	}
        
        if (write.file==T)
        {
         jpeg(file=paste(file.name,"jpg",sep="."))
         
         par(mfrow=c(3,1))
         plot((out[which(out[,1]==1),2]/1000),out[which(out[,1]==1),3], ylim=c(-10,10), type="h", xlab="position on chromosome [kb]",ylab="score", main="Chromosome 1")
         abline(h=c(0), col="black")
         abline(h=c(-2,2), col="red")
         abline(v=c(3753.687,3789.421), col="grey")
         plot((out[which(out[,1]==2),2]/1000),out[which(out[,1]==2),3], ylim=c(-10,10), type="h", xlab="position on chromosome [kb]",ylab="score", main="Chromosome 2")
         abline(h=c(0), col="black")
         abline(h=c(-2,2), col="red")
         abline(v=c(1604.786,1638.232), col="grey")
         plot((out[which(out[,1]==3),2]/1000),out[which(out[,1]==3),3], ylim=c(-10,10), type="h", xlab="position on chromosome [kb]",ylab="score", main="Chromosome 3")
         abline(h=c(0), col="black")
         abline(h=c(-2,2), col="red")
         abline(v=c(1073.708,1137.003), col="grey")
        
         dev.off()

         write.table(out,file=paste(file.name,"txt",sep="."),sep="\t",quote=F,row.names=F)         
        }
        return(out)
}


####################################################

genes.density<-function(gff, li, win=20, plot=T, write.file=F, file.name="genes.density.output", cor.method="bonferroni")

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
        if (plot==T)
        {
	 par(mfrow=c(3,1))
         plot((out[which(out[,1]==1),2]/1000),out[which(out[,1]==1),3], ylim=c(0,15), type="l", xlab="position on chromosome [kb]",ylab="-log10 of P value", main="Chromosome 1")
         abline(h=c(1.3), col="orange")
         abline(h=c(3), col="red")
         plot((out[which(out[,1]==2),2]/1000),out[which(out[,1]==2),3], ylim=c(0,15), type="l", xlab="position on chromosome [kb]",ylab="-log10 of P value", main="Chromosome 2")
         abline(h=c(1.3), col="orange")
         abline(h=c(3), col="red")
         plot((out[which(out[,1]==3),2]/1000),out[which(out[,1]==3),3], ylim=c(0,15), type="l", xlab="position on chromosome [kb]",ylab="-log10 of P value", main="Chromosome 3")
         abline(h=c(1.3), col="orange")
         abline(h=c(3), col="red")
	}
        
        if (write.file==T)
        {
         jpeg(file=paste(file.name,"jpg",sep="."))
         
         par(mfrow=c(3,1))
         plot((out[which(out[,1]==1),2]/1000),out[which(out[,1]==1),3], ylim=c(0,15), type="l", xlab="position on chromosome [kb]",ylab="-log10 of P value", main="Chromosome 1")
         abline(h=c(1.3), col="orange")
         abline(h=c(3), col="red")
         plot((out[which(out[,1]==2),2]/1000),out[which(out[,1]==2),3], ylim=c(0,15), type="l", xlab="position on chromosome [kb]",ylab="-log10 of P value", main="Chromosome 2")
         abline(h=c(1.3), col="orange")
         abline(h=c(3), col="red")
         plot((out[which(out[,1]==3),2]/1000),out[which(out[,1]==3),3], ylim=c(0,15), type="l", xlab="position on chromosome [kb]",ylab="-log10 of P value", main="Chromosome 3")
         abline(h=c(1.3), col="orange")
         abline(h=c(3), col="red")
         
         dev.off()

         write.table(out,file=paste(file.name,"txt",sep="."),sep="\t",quote=F,row.names=F)         
        }
        return(out)
}
