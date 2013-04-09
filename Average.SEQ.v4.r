plot.Average.SEQ=function(data,my.ylim=c(0.5,2),my.col=4,what="gene",my.list=NULL,para=c("cl",500,20,10,300),compute=FALSE,help=FALSE,my.log=TRUE,my.gff=gff,my.type="b",leg=TRUE,title=NULL,my.cex=1,my.lwd=1)
{
 mycol=c("purple","brown","red","orange","yellow","light green","light blue", "blue", "dark blue","black","dark grey","grey","purple","brown","pink")
 if(help==TRUE)
 {
  print("para arguments are: type,flank,nbin,fbin,dist.")
  return("finished")
 } 
 if(compute==TRUE)
 {
  data1=list()
  for(i in data)
  {
   print(i)
   data1[[i]]=Average.SEQ.all(sub_col=i,type=para[1],flank=as.numeric(para[2]),nbin=as.numeric(para[3]),fbin=as.numeric(para[4]),dist=as.numeric(para[5]),li=my.list,unlog=my.log,gffa=my.gff)
  }
 } 
 else
 {
  data1=data
 }
 
 plot(data1[[1]][,my.col],type=my.type,pch=20,col=mycol[1],ylim=my.ylim,main=title,cex=my.cex,lwd=my.lwd,ylab="A.U.",xlab="Bins")
 if(length(data1) != 1)
 {
  for(i in 2:length(data1))
  {
   lines(data1[[i]][,my.col],type=my.type,pch=20,col=mycol[i],cex=my.cex,lwd=my.lwd)
  }
 }
 if(grepl(para[1],"classic")==TRUE)
 {
  abline(v=data1[[1]][1:2,5],col="black",lwd=3)
  if(leg==TRUE)
  {
   #legend(x="bottomright",legend=(paste(data1[[1]][3,5],'bps',sep='')),cex=1.5,bty="n")
   legend(x="bottomleft",legend=(paste(data1[[1]][3,5],'bps',sep='')),cex=0.8,bty="n")
   legend(x="bottom",legend=what,cex=3,bty="n")
  }
 }
 else
 {
  abline(v=c((as.numeric(para[4])+0.5),(as.numeric(para[4])+(2*as.numeric(para[3]))+3)),col="black",lwd=3)
  abline(v=c((as.numeric(para[4])+as.numeric(para[3])+1),(as.numeric(para[4])+as.numeric(para[3])+3)),col="grey",lwd=3)
  if(leg==TRUE)
  {
   #legend(x="bottomright",legend=(paste(data1[[1]][3,5],'bps',sep='')),cex=1.5,bty="n")
   legend(x="bottomleft",legend=(paste(data1[[1]][3,5],'bps',sep='')),cex=0.8,bty="n")
  }
#  legend(x="bottom",legend=what,cex=3,bty="n")
 }
 if(compute==TRUE)
 {
  return(data1)
 }
}


test.Average.SEQ<-function(chr=1,gffa=mygff)
{
 gffg=gffa[which(gffa[,3]=="gene"),]
 gffg=gffg[which(gffg[,13]==chr),]
 print(nrow(gffg))
 if(chr==1)
 {
  out=rep(1,(2*5579133))
  chrl=5579133
 }
 if(chr==2)
 {
  out=rep(1,(2*4539804))
  chrl=4539804
 }
 if(chr==3)
 {
  out=rep(1,(2*2452883))
  chrl=2452883
 }
flat=0
if (flat!=1)
{
 for (i in 1:nrow(gffg))
 {
 
 if(gffg[i,7]=="+")
  {
   start=gffg[i,4]
   end=gffg[i,5]
   if(end > ((2*5579133)-300)){next}
   out[(start-100):(end+300)]=5
  }
  else if (gffg[i,7]=="-")
  {
   start=gffg[i,4]+chrl
   end=gffg[i,5]+chrl
   if(start < 300){next}
   out[(start-300):(end+100)]=5
  }
   
 #out[start:end]=5
 } 
}
 return(out)

}
############################
Average.SEQ.all<-function(input=list(chr1,chr2,chr3),unlog=TRUE,sub_col=NULL, gffa=mygff, flank=500, nbin=60, fbin=40, li=NULL, type=c("classic","half","set_distance"),dist=300,test_chr=NULL)
{
 norm_up=0
 norm_dwn=0
 norm_gene=0

#####FORMAT DATA INPUTS###########
 if(is.null(test_chr)==FALSE)
 {
  data=list()
  data[[1]]=input
  data[[2]]=input
  data[[3]]=input
 }
 else if(is.null(sub_col)==FALSE)
 {
  data=list()
  data[[1]]=input[[1]][,sub_col]
  data[[2]]=input[[2]][,sub_col]
  data[[3]]=input[[3]][,sub_col]
 }
 else
 {
  data=input
 }

 if(unlog==TRUE)
 {
  data[[1]]=toLINEAR(data[[1]])
  data[[2]]=toLINEAR(data[[2]])
  data[[3]]=toLINEAR(data[[3]])
 }

######FORMAT GFF##################
 gffg=gffa[which(gffa[,3]=="gene"),]
 if(is.null(test_chr)==FALSE)
 {
  gffg=gffg[which(gffg[,13] == test_chr),]
 }
 else
 {
  gffg=gffg[which(gffg[,13] %in% c(1,2,3)),]
 }
 
 if(is.null(li)==F)
 {
  gffo=gffg
  gffg=gffg[which(gffg$Name %in% li),]
  print(paste("analysing ",nrow(gffg)," genes out of ",length(li)," genes in list.",sep=''))
 }
 
#####TYPE OF PROFILE############# 

  if(grepl(type,"classic")==TRUE)
  {
   tr_out=matrix(NA,nbin,3)
  }
  else
  {
   tr_out=matrix(NA,((2*nbin)+3),3)
  }
  up_out=matrix(NA,fbin,3)
  dwn_out=matrix(NA,fbin,3)
  tr_all=data.frame(stringsAsFactors=F)

  if(is.null(flank)==FALSE)
  {
   flank_bin_size=(flank-(flank %% fbin))/fbin
   flank_bin_res=flank %% fbin
 
   flank_array_bin_size=rep(flank_bin_size,fbin)
   to_res=sample(x=fbin,size=flank_bin_res)
   flank_array_bin_size[to_res]=flank_array_bin_size[to_res]+1
  }



######################################
#MAIN LOOP
######################################

 for (i in 1:nrow(gffg))
 {

###EXTRACT COORDINATES##############
 
  length=gffg[i,5]-gffg[i,4]+1
  f1=data[[gffg[i,13]]]
  chrl=length(f1)/2

#get indexes for both orientations
#leaving start as left and end as right 
  if(gffg[i,7]=="+")
  {
   start=gffg[i,4]
   end=gffg[i,5]  
  }
  else if (gffg[i,7]=="-")
  {
   start=gffg[i,4]+chrl
   end=gffg[i,5]+chrl
  }

###EXTRACT DATA##################

  if(grepl(type,"classic")==TRUE)
  {
   tr1=f1[start:end]
   rounds=1
  }
  else if (grepl(type,"half")==TRUE)  
  {
   tr1=f1[start:(start+(length/2)-1)]
   tr2=f1[(start+(length/2)):end]
   rounds=2
  }
  else if (grepl(type,"set_distance")==TRUE)
  {
   if(length < (2*dist)){next}  
   tr1=f1[start:(start+dist-1)] 
   tr2=f1[(end-dist+1):end] 
   rounds=2
  }

#if flank == number 
#use set length
  if(is.null(flank)==FALSE)
  {    
   if (((start-flank)<0)|((end+flank)>length(f1))){next} 
   up=f1[(start-flank):(start-1)]
   dwn=f1[(end+1):(end+flank)]
  }

#if flank == null
#use intergenic region from gff
  else if (is.null(flank)==TRUE)
  {
   j=i-1
   k=i+1
   if(j==0){next}
   if(k > nrow(gffg)){next}
   if(gffg$chr[i] != gffg$chr[j]){next} 
   if(gffg$chr[i] != gffg$chr[k]){next} 
 
   if(is.null(li)==TRUE)
   {
    flank_up=(gffg$start[i]-1)-(gffg$end[j]+1)+1
    flank_dwn=(gffg$start[k]-1)-(gffg$end[i]+1)+1
   }
   else if(is.null(li)==FALSE)
   {
    flank_up=(gffo$start[which(gffo$Name==gffg$Name[i])]-1)-(gffo$end[(which(gffo$Name==gffg$Name[i])-1)]+1)+1
    flank_dwn=(gffo$start[(which(gffo$Name==gffg$Name[i])+1)]-1)-(gffo$end[which(gffo$Name==gffg$Name[i])]+1)+1
   }
   
   if (((start-flank_up)<0)|((end+flank_dwn)>length(f1))){next}
   
   up=f1[(start-flank_up):(start-1)]
   dwn=f1[(end+1):(end+flank_dwn)]
   
   if(flank_up < 0)
   {
    up=0
   }
   if(flank_dwn < 0)
   {
    dwn=0
   }
  }

#keep only genes with upstream and downstream
#intergenic regions

if((length(up)==1)|(length(dwn)==1)){next}

#head to tail for reverse strand 

 for(n in 1:rounds)
 {
  tr=get(paste("tr",n,sep=''))
  
  if((gffg[i,7]=="-"))
  {
   tr=rev(tr)
  }

#create gene bins and
#distribute residual randomply

  tr_bin_size=(length(tr)-(length(tr) %% nbin))/nbin
  tr_bin_res=length(tr) %% nbin
  array_bin_size=rep(tr_bin_size,nbin)
  to_res=sample(x=nbin,size=tr_bin_res)
  array_bin_size[to_res]=array_bin_size[to_res]+1

  l=0

#loops through bins
  if(n==1)
  {
   norm_gene=norm_gene+1
  }
  for (j in 1:nbin)
  {
   
   bin_size=array_bin_size[j]
  
   if(n==2)
   {
    m=j+nbin+3
   }
   else
   {
    m=j
   }

#set coordinates
   k=l+1
   l=k+bin_size-1

#adds up hits per base of gene i in bin j to bin j
   tr_out[m,1]=sum(tr_out[m,1],sum(tr[k:l]),na.rm=T)
#adds up number of nt of gene i in bin j to bin j 
   tr_out[m,2]=sum(tr_out[m,2],length(tr[k:l]),na.rm=T)
#adds average hits per base in bin j of gene i to bin j
   tr_out[m,3]=sum(tr_out[m,3],(sum(tr[k:l])/length(tr[k:l])),na.rm=T)
#record average hit per base in bin j of gene i
   ###tr_all[i,j]=(sum(tr[k:l])/length(tr[k:l]))
  }

  l=0
 }

#same as above for upstream and downstream regions
  l=0
  k=0
  if(gffg[i,7]=="-")
  {
   up1=rev(dwn)
   dwn1=rev(up)
   up=up1
   dwn=dwn1  
  }
  
  if(is.null(flank)==FALSE)
  {
   for (j in 1:fbin) 
   {
   
    bin_size=flank_array_bin_size[j]  

    k=l+1
    l=k+bin_size-1

    up_out[j,1]=sum(up_out[j,1],sum(up[k:l]),na.rm=T)
    up_out[j,2]=sum(up_out[j,2],length(up[k:l]),na.rm=T)
    up_out[j,3]=sum(up_out[j,3],(sum(up[k:l])/length(up[k:l])),na.rm=T)
    dwn_out[j,1]=sum(dwn_out[j,1],sum(dwn[k:l]),na.rm=T)
    dwn_out[j,2]=sum(dwn_out[j,2],length(dwn[k:l]),na.rm=T)
    dwn_out[j,3]=sum(dwn_out[j,3],(sum(dwn[k:l])/length(dwn[k:l])),na.rm=T)
   }
  }
  
  else if (is.null(flank)==TRUE) 
  {
   if(length(up) != 1)
   {
    up_bin_size=(length(up)-(length(up) %% fbin))/fbin
    up_bin_res=length(up) %% fbin
    up_array_bin_size=rep(up_bin_size,fbin)
    to_res=sample(x=fbin,size=up_bin_res)
    up_array_bin_size[to_res]=up_array_bin_size[to_res]+1
   
    lu=0
    ku=0
  
    for (j in 1:fbin) 
    {
    
     up_bin_size=up_array_bin_size[j]  

     ku=lu+1
     lu=ku+up_bin_size-1

     up_out[j,1]=sum(up_out[j,1],sum(up[ku:lu]),na.rm=T)
     up_out[j,2]=sum(up_out[j,2],length(up[ku:lu]),na.rm=T)
     up_out[j,3]=sum(up_out[j,3],(sum(up[ku:lu])/length(up[ku:lu])),na.rm=T)
    }
   }
   else
   {
    norm_up=norm_up+1
   }
   if(length(dwn) != 1)
   {
    dwn_bin_size=(length(dwn)-(length(dwn) %% fbin))/fbin
    dwn_bin_res=length(dwn) %% fbin
    dwn_array_bin_size=rep(dwn_bin_size,fbin)
    to_res=sample(x=fbin,size=dwn_bin_res)
    dwn_array_bin_size[to_res]=dwn_array_bin_size[to_res]+1

    ld=0
    kd=0

    for (j in 1:fbin) 
    {
    
     dwn_bin_size=dwn_array_bin_size[j]  

     kd=ld+1
     ld=kd+dwn_bin_size-1

     dwn_out[j,1]=sum(dwn_out[j,1],sum(dwn[kd:ld]),na.rm=T)
     dwn_out[j,2]=sum(dwn_out[j,2],length(dwn[kd:ld]),na.rm=T)
     dwn_out[j,3]=sum(dwn_out[j,3],(sum(dwn[kd:ld])/length(dwn[kd:ld])),na.rm=T)
    }
   }
   else
   {
    norm_dwn=norm_dwn+1
   }
  }
 }

 ####creates output
 out=rbind(up_out,tr_out,dwn_out)
 out=cbind(out,round((out[,3]/out[1,3]),2))
 bars=rep(NA,nrow(out))
 bars[1]=fbin+0.5
 bars[2]=fbin+nbin+0.5
 if(is.null(flank)==FALSE)
 {
  bars[3]=flank
 }
 else
 {
  bars[3]=0
 }
 out=cbind(out,bars)
#print(out)
 out=cbind(out,round((out[,3]/out[fbin,3]),2)) 
#return(list(out,tr_all))
 print(norm_up)
 print(norm_dwn)
 print(norm_gene)
 return(out)
}








