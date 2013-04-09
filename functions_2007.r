####################################################################################
viewSEQ.dots<-function(c1=chr1, c2=chr2, c3=chr3, range=c(1478,1490),scale=c(-4,20), highlights=c(1487000,1489000,1481000,1485000), gff=mygff,chr=1, sep=F,...){
grid.newpage()
if (chr==1){n=c1}
else if (chr==2){n=c2}
else if (chr==3){n=c3}
highlight1=list()
highlight1$coord=as.integer(highlights[1])
highlight1$strand="+"

highlight2=list()
highlight2$coord=as.integer(highlights[2])
highlight2$strand="+"

highlight3=list()
highlight3$coord=as.integer(highlights[3])
highlight3$strand="-"

highlight4=list()
highlight4$coord=as.integer(highlights[4])
highlight4$strand="-"

#str(paste("chr",chr,sep=""))
MyplotAlong(highlight1=highlight1,highlight2=highlight2, highlight3=highlight3, highlight4=highlight4, y = as.matrix(n), chr = chr, coord = range * 1000, ylim=scale, gff = mygff, probeAnno=probeAnno, what=c("dots"),sepPlots=sep, featureExclude = c("LTR","misc_RNA","new","new_antisense", "gap","3'UTR", "5'UTR", "CDS_motif", "gap", "mRNA", "polyA_site", "promoter", "rep_origin", "repeat_unit","polyA_signal","BLASTN_HIT", "trans"))

gc()
}
####################################################################################
viewSEQ.sep<-function(c1=chr1, c2=chr2, c3=chr3, range=c(1478,1490),scale=c(-4,20), highlights=c(1487000,1489000,1481000,1485000), gff=mygff,chr=1, sep=T,type="dots",...){
grid.newpage()
if (chr==1){n=c1}
else if (chr==2){n=c2}
else if (chr==3){n=c3}
highlight1=list()
highlight1$coord=as.integer(highlights[1])
highlight1$strand="+"

highlight2=list()
highlight2$coord=as.integer(highlights[2])
highlight2$strand="+"

highlight3=list()
highlight3$coord=as.integer(highlights[3])
highlight3$strand="-"

highlight4=list()
highlight4$coord=as.integer(highlights[4])
highlight4$strand="-"

#str(paste("chr",chr,sep=""))
MyplotAlong(highlight1=highlight1,highlight2=highlight2, highlight3=highlight3, highlight4=highlight4, y = as.matrix(n), chr = chr, coord = range * 1000, gff = mygff, probeAnno=probeAnno, what=type,sepPlots=sep, featureExclude = c("LTR","misc_RNA","new","new_antisense", "gap","3'UTR", "5'UTR", "CDS_motif", "gap", "mRNA", "polyA_site", "promoter","rep_origin", "repeat_unit","polyA_signal","BLASTN_HIT", "trans","utr_intron"))

gc()
}
###############################################################################################
viewSEQ.heat<-function(c1=chr1, c2=chr2, c3=chr3, range=c(1478,1490),scale=c(-4,20), highlights=c(1487000,1489000,1481000,1485000), gff=mygff,chr=1, sep=F,...){
grid.newpage()
if (chr==1){n=c1}
else if (chr==2){n=c2}
else if (chr==3){n=c3}
highlight1=list()
highlight1$coord=as.integer(highlights[1])
highlight1$strand="+"

highlight2=list()
highlight2$coord=as.integer(highlights[2])
highlight2$strand="+"

highlight3=list()
highlight3$coord=as.integer(highlights[3])
highlight3$strand="-"

highlight4=list()
highlight4$coord=as.integer(highlights[4])
highlight4$strand="-"

#str(paste("chr",chr,sep=""))
MyplotAlong(highlight1=highlight1,highlight2=highlight2, highlight3=highlight3, highlight4=highlight4, y = as.matrix(n), chr = chr, coord = range * 1000, ylim = scale, gff = mygff, probeAnno=probeAnno, what="heatmap",sepPlots=sep, featureExclude = c("LTR","misc_RNA","new","new_antisense", "gap","3'UTR", "5'UTR", "CDS_motif", "gap", "mRNA", "polyA_site", "promoter","rep_origin", "repeat_unit","polyA_signal","BLASTN_HIT", "trans","utr_intron"))

gc()
}
####################################################################################
viewGEN.sep<-function(gene='SPAC1250.05',zoom=0, c1=chr1, c2=chr2, c3=chr3,flank=2,scale=c(-4,20), gff=mygff, sep=T,high=c(0,0),col=0,verbose=T,what="dots",...){
grid.newpage()
t=getc(gene,all=T)
chr=as.numeric(t[2,2])
if (chr==1){n=c1}
else if (chr==2){n=c2}
else if (chr==3){n=c3}

if (col != 0)
{
 n=n[,col]
}

start=floor(as.numeric(t[2,4])/1000)-flank
end=ceiling(as.numeric(t[2,5])/1000)+flank

if (length(zoom)==1)
{
range=c(start,end)
}
else
{
range=zoom
}

if(length(high)==4)
{
 highlight1=list()
 highlight1$coord=high[1]
 highlight1$strand="+"

 highlight2=list()
 highlight2$coord=high[2]
 highlight2$strand="+"

 highlight3=list()
 highlight3$coord=high[3]
 highlight3$strand="-"

 highlight4=list()
 highlight4$coord=high[4]
 highlight4$strand="-"
}

else
{
 highlight1=list()
 highlight1$coord=as.numeric(t[2,4])
 highlight1$strand="+"

 highlight2=list()
 highlight2$coord=as.numeric(t[2,5])
 highlight2$strand="+"

 highlight3=list()
 highlight3$coord=as.numeric(t[1,4])
 highlight3$strand="-"

 highlight4=list()
 highlight4$coord=as.numeric(t[1,5])
 highlight4$strand="-"
}

if (verbose==T)
{
print(chr)
print(start)
print(end)
}
#str(paste("chr",chr,sep=""))
MyplotAlong(highlight1=highlight1, highlight2=highlight2, highlight3=highlight3, highlight4=highlight4, y = as.matrix(n), chr = chr, coord = range * 1000, gff = mygff, probeAnno=probeAnno, what=what,sepPlots=sep, featureExclude = c("LTR","misc_RNA","new","new_antisense", "gap","3'UTR", "5'UTR", "CDS_motif", "gap", "mRNA", "polyA_site", "promoter","rep_origin", "repeat_unit","polyA_signal","BLASTN_HIT", "trans","utr_intron"))
gc(verbose=F)
}
###############################################################################################
#swapChr=function(chr1=chr1,chr2=chr2,chr3=chr3,chr1b=chr1b,chr2b=chr2b,chr3b=chr3b)
#{
#chr1a=chr1b
#chr1b=chr1
#chr1=chr1a
#chr2a=chr2b
#chr2b=chr2
#chr2=chr2a
#chr3a=chr3b
#chr3b=chr3
#chr3=chr3a
#return(chr1,chr2,chr3,chr1b,chr2b,chr3b)
#}
#swapPrep=function(chr1b=chr1,chr2b=chr2,chr3b=chr3)
#{
#return(chr1b,chr2b,chr3b)
#}
#####################################################################################
set.mygff<-function(){
mygff = which(gff[, "start"] <= gff[, "end"])
mygff=gff[mygff,]
return(mygff)
}
################################################################################################
getc.trans=read.delim("allNames_030811.TXT2", stringsAsFactors=F)
################################################################################################
getc=function(gene,all=F)
{

gene1=getc.trans[grep(gene, getc.trans[,2]),1]
if(length(gene1) > 1)
{
gene1=getc.trans[which(getc.trans[,2]==gene),1]
}

if (all==T)
{
out=mygff[which((mygff[,"Name"]==gene1)),c("Name","chr","feature","start","end","strand","attributes")]
}
else if (all==F)
{
out=mygff[which((mygff[,"feature"]=="gene") & (mygff[,"Name"]==gene1)),c("Name","chr","start","end","strand","attributes")]
}
return(out)
}
################################################################################################
####MODIFIED TILINGARRAY FUNCTION###############################################################
################################################################################################
MyplotAlongSeq<-
function (segObj, y, probeAnno, gff, isDirectHybe = FALSE, what = c("dots","heatmap"), 
    chr, coord, highlight1,highlight2, colors, doLegend = FALSE, featureExclude = c("chromosome", 
        "nucleotide_match", "insertion"), featureColorScheme = 1, 
    extras, rowNamesHeatmap, rowNamesExtras, ylab, ylabExtras, 
    main, colHeatmap = colorRamp(brewer.pal(9, "YlOrBr")), colExtras = colorRamp(brewer.pal(9, 
        "Reds")), sepPlots = FALSE, reOrder = TRUE, mixScale=FALSE, scaleD=NULL, scaleU=NULL, ...) 
{
    VP = c(title = 0.4, "expr+" = 5, "gff+" = 1, coord = 1, "gff-" = 1, 
        "expr-" = 5, legend = 0.4)
    if (sepPlots) {
        if (!missing(y)) 
            n <- ncol(y)
        if (!missing(segObj)) {
            if (is.environment(segObj)) {
                segmentationObjectName = paste(chr, "+", sep = ".")
                if (segmentationObjectName %in% ls(segObj)) {
                  s <- get(segmentationObjectName, segObj)
                  n <- ncol(s@y)
                }
                else {
                  dat = get(paste(chr, strand, "dat", sep = "."), 
                    segObj)
                  n <- ncol(dat$y)
                }
            }
        }
        if (reOrder) 
            ordering = seq(n, 1, by = -1)
        else ordering = seq(1:n)
        if (n < 4) {
            exprw <- exprc <- NULL
            for (i in 1:n) {
                exprw <- c(exprw, paste("expr", i, "+", sep = ""))
                exprc <- c(exprc, paste("expr", i, "-", sep = ""))
            }
            VPnames <- c("title", exprw, "gff+", "coord", "gff-", 
                exprc, "legend")
            VP = c(0.8, rep(5, n), 1, 1, 1, rep(5, n), 0.4)
            names(VP) <- VPnames
        }
        else {
            cat("More than 4 arrays, plotting averages.\n")
            sepPlots = FALSE
        }
    }
    if (!missing(extras)) {
        indgff <- grep("gff\\+", names(VP))
        indlegend <- grep("legend", names(VP))
        VP <- c(VP[1:(indgff - 1)], "extras+" = 1, VP[indgff:(indlegend - 
            1)], "extras-" = 1, VP[indlegend:length(VP)])
    }
    if (!doLegend) 
        VP = VP[-which(names(VP) == "legend")]
    if (missing(gff)) 
        VP = VP[-which(names(VP) == "gff+" | names(VP) == "gff-")]
    defaultColors = c("+" = "#00441b", "-" = "#081d58", duplicated = "red", 
        cp = "#555555", ci = "#777777", highlight = "red", threshold = "grey")
    if (!missing(colors)) {
        mt = match(names(colors), names(defaultColors))
        if (any(is.na(mt))) 
            stop(paste("Cannot use color specification for", 
                names(colors)[is.na(mt)]))
        defaultColors[mt] = colors
    }
    colors = defaultColors
    if (!missing(y)) {
        if (missing(probeAnno)) 
            stop("If 'y' is specified, 'probeAnno' must also be specified.")
        if (!missing(segObj)) 
            stop("If 'y' is specified, 'segObj' must not be specified.")
    }
    else {
        if (missing(segObj)) 
            stop("Please specify either 'y' or 'segObj'")
    }
    pushViewport(viewport(width = 0.85, height = 0.95))
    pushViewport(viewport(layout = grid.layout(length(VP), 1, 
        height = VP)))
    for (i in 1:2) {
    #for (i in 1) {
        strand = c("+", "-")[i]
        threshold = as.numeric(NA)
        if (!missing(y)) {
            stopifnot(is.matrix(y))
            index = get(paste(chr, strand, "index", sep = "."), 
                envir = probeAnno)
            sta = get(paste(chr, strand, "start", sep = "."), 
                envir = probeAnno)
            end = get(paste(chr, strand, "end", sep = "."), envir = probeAnno)
            if (!missing(extras)) 
                dat = list(x = (sta + end)/2, y = y[index, , 
                  drop = FALSE], flag = get(paste(chr, strand, 
                  "unique", sep = "."), envir = probeAnno), extras = extras[index, 
                  , drop = FALSE])
            else dat = list(x = (sta + end)/2, y = y[index, , 
                drop = FALSE], flag = get(paste(chr, strand, 
                "unique", sep = "."), envir = probeAnno))
            stopifnot(is.numeric(dat$flag))
            lengthChr = end[length(end)]
        }
        else {
            if (!is.environment(segObj)) 
                stop("'segObj' must be an environment.")
            segmentationObjectName = paste(chr, strand, sep = ".")
            if (segmentationObjectName %in% ls(segObj)) {
                s = get(segmentationObjectName, segObj)
                if (!inherits(s, "segmentation")) 
                  stop(sprintf("'%s' must be of class'segmentation'.", 
                    segmentationObjectName))
                if (is.na(s@nrSegments)) 
                  stop(sprintf("Slot 'nrSegments' of '%s' must not be NA.", 
                    segmentationObjectName))
                bp = s@breakpoints[[s@nrSegments]]
                dat = list(x = s@x, y = s@y, flag = s@flag, estimate = bp[, 
                  "estimate"])
                if ("upper" %in% colnames(bp)) 
                  dat$upper = bp[, "upper"]
                if ("lower" %in% colnames(bp)) 
                  dat$lower = bp[, "lower"]
            }
            else {
                dat = get(paste(chr, strand, "dat", sep = "."), 
                  segObj)
                stopifnot(all(c("start", "end", "unique", "ss") %in% 
                  names(dat)))
                dat$x = (dat$start + dat$end)/2
                dat$flag = dat$unique
                lengthChr = dat$end[length(dat$end)]
                if ("segScore" %in% ls(segObj)) {
                  sgs = get("segScore", segObj)
                  sgs = sgs[sgs$chr == chr & sgs$strand == strand, 
                    c("start", "end")]
                }
                else {
                  stop("This option is deprecated")
                }
                dat$estimate = (sgs$start[-1] + sgs$end[-length(sgs$end)])/2
                if ("theThreshold" %in% ls(segObj)) 
                  threshold = get("theThreshold", segObj)
            }
        }
        if (missing(coord)) 
            coord = c(1, lengthChr)
        vpr = which(names(VP) == sprintf("expr%s", strand))
        switch(what, dots = {
############################           

 		if (sepPlots) {
                ylimdata = quantile(as.vector(dat[["y"]][dat[["x"]] >= 
                  coord[1] & dat[["x"]] <= coord[2], ]), c(0, 
                  1), na.rm = TRUE)
                ylim = ylimdata
                #ylim=ylim
                if (missing(ylab)) 
                  ylab = colnames(dat$y)
                if (length(ylab) == 1) 
                  ylab = rep(ylab, n)
                
                for (j in seq(1:n))
                {
                 datj <- dat
                 datj$y <- dat$y[, ordering[j]]
                 if (missing(ylab)) 
                   ylab = colnames(dat$y)
                 vpr = which(names(VP) == sprintf(paste("expr", j, "%s", sep = ""), strand))

                  #plotSegmentationDots(datj, xlim = coord, ylim = ylim, 
                 MyplotSegmentationDots(datj, xlim = coord, ylim = ylim,ylab = ylab[ordering[j]], chr = chr, strand = ifelse(isDirectHybe, 
                 otherStrand(strand), strand), vpr = vpr, colors = colors, sepPlots = sepPlots,  ...)
                }
            }
            #else plotSegmentationDots(dat, xlim = coord, ylab = ylab, 
              #else MyplotSegmentationDots(dat, xlim = coord, ylab = ylab,
               else if (i==1) MyplotSegmentationDots(dat, xlim = coord, ylab = ylab,
                chr = chr, strand = ifelse(isDirectHybe, otherStrand(strand), 
                  strand), vpr = vpr, colors = colors, sepPlots = sepPlots, 
                ...)
#added sam
            if (!missing(extras) & !missing(y)) {
                vpr2 = which(names(VP) == sprintf("extras%s", 
                  strand))
                dat$y = dat$extras[, , drop = FALSE]
                #plotSegmentationDots(dat, xlim = coord, chr = chr, 
                MyplotSegmentationDots(dat, xlim = coord, chr = chr, 
                  strand = ifelse(isDirectHybe, otherStrand(strand), 
                    strand), vpr = vpr2, colors = colors, colHeatmap = colExtras, 
                  ylab = ylabExtras, rowNames = rowNamesExtras, 
                  ...)
            }
        }, heatmap = {
            #plotSegmentationHeatmap(dat, xlim = coord, rowNames = rowNamesHeatmap,
        if (i==1){    plotSegmentationHeatmap(dat, xlim = coord, rowNames = rowNamesHeatmap, 
                chr = chr, strand = ifelse(isDirectHybe, otherStrand(strand), 
                  strand), vpr = vpr, colors = colors, ylab = ylab, 
                colHeatmap = colHeatmap, ...)
            }
#added by Sam
            if (!missing(extras) & !missing(y)) {
                vpr2 = which(names(VP) == sprintf("extras%s", 
                  strand))
                dat$y = dat$extras[, , drop = FALSE]
                plotSegmentationHeatmap(dat, xlim = coord, chr = chr, 
                  strand = ifelse(isDirectHybe, otherStrand(strand), 
                    strand), vpr = vpr2, colors = colors, colHeatmap = colExtras, 
                  ylab = ylabExtras, rowNames = rowNamesExtras, 
                  ...)
            }
        }, stop(sprintf("Invalid value '%s' for argument 'what'", 
            what)))
        if (!missing(gff)) 
            plotFeatures(gff = gff, chr = chr, xlim = coord, 
                strand = strand, featureExclude = featureExclude, 
                featureColorScheme = featureColorScheme, vpr = which(names(VP) == 
                  sprintf("gff%s", strand)), ...)
    }
    pushViewport(dataViewport(xData = coord, yscale = c(-0.4, 
        0.8), extension = 0, layout.pos.col = 1, layout.pos.row = which(names(VP) == 
        "coord")))
    grid.lines(coord, c(0, 0), default.units = "native")
    tck = tilingArray:::alongChromTicks(coord)
    grid.text(label = formatC(tck, format = "d"), x = tck, y = 0.2, 
        just = c("centre", "bottom"), gp = gpar(cex = 0.6), default.units = "native")
    grid.segments(x0 = tck, x1 = tck, y0 = -0.17, y1 = 0.17, 
        default.units = "native")
##############
   if (!missing(highlight1)) {
        mt = (match(highlight1$strand, c("-", "+")) - 1.5) * 2
        co = highlight1$coord
        if (is.na(mt) || !is.numeric(co)) 
            stop("Invalid parameter 'highlight1'.")
        grid.segments(x0 = co, x1 = co, y0 = c(0, 0), y1 = c(0.4, 
            0.4) * mt, default.units = "native", gp = gpar(col = colors["highlight"], 
            lwd = 2))
    }
##############
if (!missing(highlight2)) {
        mt = (match(highlight2$strand, c("-", "+")) - 1.5) * 2
        co = highlight2$coord
        if (is.na(mt) || !is.numeric(co)) 
            stop("Invalid parameter 'highlight2'.")
        grid.segments(x0 = co, x1 = co, y0 = c(0, 0), y1 = c(0.4, 
            0.4) * mt, default.units = "native", gp = gpar(col = colors["highlight"], 
            lwd = 2))
    }
##############
    popViewport()
    pushViewport(viewport(layout.pos.col = 1, layout.pos.row = which(names(VP) == 
        "title")))
    grid.text(label = paste("Chr ", chr, sep = ""), x = 0.5, 
        y = 1, just = "centre", gp = gpar(cex = 1))
    if (!missing(main)) 
        grid.text(label = main, x = 0.05, y = 1, just = "centre", 
            gp = gpar(cex = 1))
    popViewport()
    if (doLegend) 
        plotAlongChromLegend(which(names(VP) == "legend"), featureColorScheme = featureColorScheme, 
            featureExclude = featureExclude)
    popViewport(2)
}

################################################################################
################################################################################

MyplotSegmentationDots <- function (dat, xlim, ylim, ylab, threshold = NA, chr = 1, strand = "+", 
    vpr, colors, main, pointSize = unit(0.6, "mm"), showConfidenceIntervals = TRUE, 
    sepPlots = FALSE, cexAxisLabel = 1, cexAxis = 1, ...) 
{
    endVP = FALSE
    if (missing(vpr)) {
        endVP = TRUE
        vpr = newVP(main = main, dataPanelHeight = 1, vpHeight = 0.95, 
            titleOffSet = -0.9)
    }
    if (is.matrix(dat$y) & !sepPlots) 
        dat$y = rowMeans(dat$y)
    stopifnot(length(dat$y) == length(dat$x), length(dat$flag) == 
        length(dat$x))
    xorg = dat$x
    if (missing(xlim)) {
        xlim = range(dat$x, na.rm = TRUE)
    }
    else {
        sel = (dat$x >= xlim[1]) & (dat$x <= xlim[2])
        dat$x = dat$x[sel]
        dat$y = dat$y[sel]
        dat$flag = dat$flag[sel]
    }
    if (!is.na(threshold)) {
        dat$y = dat$y - threshold
        if (!missing(ylim)) 
            ylim = ylim - threshold
    }
    if (missing(ylim)) 
        ylim = quantile(dat$y, c(0, 1), na.rm = TRUE)
    if (!missing(ylab)) {
        pushViewport(dataViewport(xData = xlim, yData = ylim, 
            extension = 0, clip = "off", layout.pos.col = 1, 
            layout.pos.row = vpr))
        grid.yaxis(gp = gpar(cex = cexAxis), ...)#removed(added .log)
        grid.text(ylab, x = -0.075, y = 0.5, rot = 90, gp = gpar(cex = cexAxisLabel), 
            ...)
        pushViewport(dataViewport(xData = xlim, yData = ylim, 
            extension = 0, clip = "on", layout.pos.col = 1, layout.pos.row = vpr))
    }
    else {
        pushViewport(dataViewport(xData = xlim, yData = ylim, 
            extension = 0, clip = "off", layout.pos.col = 1, 
            layout.pos.row = vpr))
        grid.yaxis(gp = gpar(cex = cexAxis))#removed(added.log)
        pushViewport(dataViewport(xData = xlim, yData = ylim, 
            extension = 0, clip = "on", layout.pos.col = 1, layout.pos.row = vpr))
    }
    defaultColors = c("+" = "#00441b", "-" = "#081d58", duplicated = "grey",
       cp = "#555555", ci = "#777777", highlight = "red", threshold = "grey")
    if (!missing(colors)) {
        mt = match(names(colors), names(defaultColors))
        if (any(is.na(mt))) 
            stop(paste("Cannot use color specification for", 
                names(colors)[is.na(mt)]))
        defaultColors[mt] = colors
    }
    colors = defaultColors
    ord = c(which(dat$flag != 0), which(dat$flag == 0))
    colo = ifelse(dat$flag[ord] == 0, colors[strand], colors["duplicated"])
    if (!is.na(threshold)) 
        grid.lines(y = unit(0, "native"), gp = gpar(col = colors["threshold"]))
    sel = ((xorg[dat$estimate] >= xlim[1]) & (xorg[dat$estimate] <= 
        xlim[2]))
    mySeg = function(j, what) grid.segments(x0 = unit(j, "native"), 
        x1 = unit(j, "native"), y0 = unit(0.1, "npc"), y1 = unit(0.9, 
            "npc"), gp = gpar(col = colors[what], lty = c(cp = 1, 
            ci = 2)[what]))
    if (!is.null(dat$estimate) & sum(sel) > 0) 
        mySeg(xorg[dat$estimate][sel], "cp")
    if (showConfidenceIntervals & sum(sel) > 0) {
        if (!is.null(dat$upper)) 
            mySeg(xorg[dat$upper][sel], "ci")
        if (!is.null(dat$lower)) 
            mySeg(xorg[dat$lower][sel], "ci")
    }
    grid.points(dat$x[ord], dat$y[ord], pch = 20, size = pointSize, 
        gp = gpar(col = colo))
    popViewport(2)
    if (endVP) 
        popViewport(2)
}
#################################################################################################
MyplotAlong<-
function (segObj, y, probeAnno, gff, isDirectHybe = FALSE, what = c("dots","heatmap"), 
    chr, coord, highlight1,highlight2,highlight3,highlight4, colors, doLegend = FALSE, featureExclude = c("chromosome", 
        "nucleotide_match", "insertion"), featureColorScheme = 1, 
    extras, rowNamesHeatmap, rowNamesExtras, ylab, ylabExtras, 
    main, colHeatmap = colorRamp(brewer.pal(9, "YlGnBu")), colExtras = colorRamp(brewer.pal(9, 
        "Reds")), sepPlots = FALSE, reOrder = TRUE, ...) 
{
    VP = c(title = 0.4, "expr+" = 5, "gff+" = 1, coord = 1, "gff-" = 1, 
        "expr-" = 5, legend = 0.4)
    if (sepPlots) {
        if (!missing(y)) 
            n <- ncol(y)
        if (!missing(segObj)) {
            if (is.environment(segObj)) {
                segmentationObjectName = paste(chr, "+", sep = ".")
                if (segmentationObjectName %in% ls(segObj)) {
                  s <- get(segmentationObjectName, segObj)
                  n <- ncol(s@y)
                }
                else {
                  dat = get(paste(chr, strand, "dat", sep = "."), 
                    segObj)
                  n <- ncol(dat$y)
                }
            }
        }
        if (reOrder) 
            ordering = seq(n, 1, by = -1)
        else ordering = seq(1:n)
        if (n < 4) {
            exprw <- exprc <- NULL
            for (i in 1:n) {
                exprw <- c(exprw, paste("expr", i, "+", sep = ""))
                exprc <- c(exprc, paste("expr", i, "-", sep = ""))
            }
            VPnames <- c("title", exprw, "gff+", "coord", "gff-", 
                exprc, "legend")
            VP = c(0.8, rep(5, n), 1, 1, 1, rep(5, n), 0.4)
            names(VP) <- VPnames
        }
        else {
            cat("More than 3 arrays, plotting averages.\n")
            sepPlots = FALSE
        }
    }
    if (!missing(extras)) {
        indgff <- grep("gff\\+", names(VP))
        indlegend <- grep("legend", names(VP))
        VP <- c(VP[1:(indgff - 1)], "extras+" = 1, VP[indgff:(indlegend - 
            1)], "extras-" = 1, VP[indlegend:length(VP)])
    }
    if (!doLegend) 
        VP = VP[-which(names(VP) == "legend")]
    if (missing(gff)) 
        VP = VP[-which(names(VP) == "gff+" | names(VP) == "gff-")]
    defaultColors = c("+" = "#00441b", "-" = "#081d58", duplicated = "red", 
        cp = "#555555", ci = "#777777", highlight = "red", threshold = "grey")
    if (!missing(colors)) {
        mt = match(names(colors), names(defaultColors))
        if (any(is.na(mt))) 
            stop(paste("Cannot use color specification for", 
                names(colors)[is.na(mt)]))
        defaultColors[mt] = colors
    }
    colors = defaultColors
    if (!missing(y)) {
        if (missing(probeAnno)) 
            stop("If 'y' is specified, 'probeAnno' must also be specified.")
        if (!missing(segObj)) 
            stop("If 'y' is specified, 'segObj' must not be specified.")
    }
    else {
        if (missing(segObj)) 
            stop("Please specify either 'y' or 'segObj'")
    }
    pushViewport(viewport(width = 0.85, height = 0.95))
    pushViewport(viewport(layout = grid.layout(length(VP), 1, 
        height = VP)))
    for (i in 1:2) {
    #for (i in 1) {
        strand = c("+", "-")[i]
        threshold = as.numeric(NA)
        if (!missing(y)) {
            stopifnot(is.matrix(y))
            index = get(paste(chr, strand, "index", sep = "."), 
                envir = probeAnno)
            sta = get(paste(chr, strand, "start", sep = "."), 
                envir = probeAnno)
            end = get(paste(chr, strand, "end", sep = "."), envir = probeAnno)
            if (!missing(extras)) 
                dat = list(x = (sta + end)/2, y = y[index, , 
                  drop = FALSE], flag = get(paste(chr, strand, 
                  "unique", sep = "."), envir = probeAnno), extras = extras[index, 
                  , drop = FALSE])
            else dat = list(x = (sta + end)/2, y = y[index, , 
                drop = FALSE], flag = get(paste(chr, strand, 
                "unique", sep = "."), envir = probeAnno))
            stopifnot(is.numeric(dat$flag))
            lengthChr = end[length(end)]
        }
        else {
            if (!is.environment(segObj)) 
                stop("'segObj' must be an environment.")
            segmentationObjectName = paste(chr, strand, sep = ".")
            if (segmentationObjectName %in% ls(segObj)) {
                s = get(segmentationObjectName, segObj)
                if (!inherits(s, "segmentation")) 
                  stop(sprintf("'%s' must be of class'segmentation'.", 
                    segmentationObjectName))
                if (is.na(s@nrSegments)) 
                  stop(sprintf("Slot 'nrSegments' of '%s' must not be NA.", 
                    segmentationObjectName))
                bp = s@breakpoints[[s@nrSegments]]
                dat = list(x = s@x, y = s@y, flag = s@flag, estimate = bp[, 
                  "estimate"])
                if ("upper" %in% colnames(bp)) 
                  dat$upper = bp[, "upper"]
                if ("lower" %in% colnames(bp)) 
                  dat$lower = bp[, "lower"]
            }
            else {
                dat = get(paste(chr, strand, "dat", sep = "."), 
                  segObj)
                stopifnot(all(c("start", "end", "unique", "ss") %in% 
                  names(dat)))
                dat$x = (dat$start + dat$end)/2
                dat$flag = dat$unique
                lengthChr = dat$end[length(dat$end)]
                if ("segScore" %in% ls(segObj)) {
                  sgs = get("segScore", segObj)
                  sgs = sgs[sgs$chr == chr & sgs$strand == strand, 
                    c("start", "end")]
                }
                else {
                  stop("This option is deprecated")
                }
                dat$estimate = (sgs$start[-1] + sgs$end[-length(sgs$end)])/2
                if ("theThreshold" %in% ls(segObj)) 
                  threshold = get("theThreshold", segObj)
            }
        }
        if (missing(coord)) 
            coord = c(1, lengthChr)
        vpr = which(names(VP) == sprintf("expr%s", strand))
        switch(what, dots = {
            if (sepPlots) {
                ylimdata = quantile(as.vector(dat[["y"]][dat[["x"]] >= 
                  coord[1] & dat[["x"]] <= coord[2], ]), c(0, 
                  1), na.rm = TRUE)
                ylim = ylimdata
                #ylim=ylim
                if (missing(ylab)) 
                  ylab = colnames(dat$y)
                if (length(ylab) == 1) 
                  ylab = rep(ylab, n)
                for (j in seq(1:n)) {
                  datj <- dat
                  datj$y <- dat$y[, ordering[j]]
                  if (missing(ylab)) 
                    ylab = colnames(dat$y)
                  vpr = which(names(VP) == sprintf(paste("expr", 
                    j, "%s", sep = ""), strand))
                  #plotSegmentationDots(datj, xlim = coord, ylim = ylim, 
                    MyplotSegmentationDots(datj, xlim = coord, ylim = ylim,
                    ylab = ylab[ordering[j]], chr = chr, strand = ifelse(isDirectHybe, 
                      otherStrand(strand), strand), vpr = vpr, 
                    colors = colors, sepPlots = sepPlots,  ...)
                }
            }
            #else plotSegmentationDots(dat, xlim = coord, ylab = ylab, 
              else MyplotSegmentationDots(dat, xlim = coord, ylab = ylab,
               #else if (i==1) MyplotSegmentationDots(dat, xlim = coord, ylab = ylab,#added sam
                chr = chr, strand = ifelse(isDirectHybe, otherStrand(strand), 
                  strand), vpr = vpr, colors = colors, sepPlots = sepPlots, 
                ...)
            if (!missing(extras) & !missing(y)) {
                vpr2 = which(names(VP) == sprintf("extras%s", 
                  strand))
                dat$y = dat$extras[, , drop = FALSE]
                #plotSegmentationDots(dat, xlim = coord, chr = chr, 
                MyplotSegmentationDots(dat, xlim = coord, chr = chr, 
                  strand = ifelse(isDirectHybe, otherStrand(strand), 
                    strand), vpr = vpr2, colors = colors, colHeatmap = colExtras, 
                  ylab = ylabExtras, rowNames = rowNamesExtras, 
                  ...)
            }
        }, heatmap = {
            plotSegmentationHeatmap(dat, xlim = coord, rowNames = rowNamesHeatmap, 
                chr = chr, strand = ifelse(isDirectHybe, otherStrand(strand), 
                  strand), vpr = vpr, colors = colors, ylab = ylab, 
                colHeatmap = colHeatmap, ...)
            if (!missing(extras) & !missing(y)) {
                vpr2 = which(names(VP) == sprintf("extras%s", 
                  strand))
                dat$y = dat$extras[, , drop = FALSE]
                plotSegmentationHeatmap(dat, xlim = coord, chr = chr, 
                  strand = ifelse(isDirectHybe, otherStrand(strand), 
                    strand), vpr = vpr2, colors = colors, colHeatmap = colExtras, 
                  ylab = ylabExtras, rowNames = rowNamesExtras, 
                  ...)
            }
        }, stop(sprintf("Invalid value '%s' for argument 'what'", 
            what)))
        if (!missing(gff)) 
            plotFeatures(gff = gff, chr = chr, xlim = coord, 
                strand = strand, featureExclude = featureExclude, 
                featureColorScheme = featureColorScheme, vpr = which(names(VP) == 
                  sprintf("gff%s", strand)), ...)
    }
    pushViewport(dataViewport(xData = coord, yscale = c(-0.4, 
        0.8), extension = 0, layout.pos.col = 1, layout.pos.row = which(names(VP) == 
        "coord")))
    grid.lines(coord, c(0, 0), default.units = "native")
    tck = tilingArray:::alongChromTicks(coord)
    grid.text(label = formatC(tck, format = "d"), x = tck, y = 0.2, 
        just = c("centre", "bottom"), gp = gpar(cex = 0.6), default.units = "native")
    grid.segments(x0 = tck, x1 = tck, y0 = -0.17, y1 = 0.17, 
        default.units = "native")
##############
   if (!missing(highlight1)) {
        mt = (match(highlight1$strand, c("-", "+")) - 1.5) * 2
        co = highlight1$coord
        if (is.na(mt) || !is.numeric(co)) 
            stop("Invalid parameter 'highlight1'.")
        grid.segments(x0 = co, x1 = co, y0 = c(0, 0), y1 = c(0.4, 
            0.4) * mt, default.units = "native", gp = gpar(col = colors["highlight"], 
            lwd = 2))
    }
##############
if (!missing(highlight2)) {
        mt = (match(highlight2$strand, c("-", "+")) - 1.5) * 2
        co = highlight2$coord
        if (is.na(mt) || !is.numeric(co)) 
            stop("Invalid parameter 'highlight2'.")
        grid.segments(x0 = co, x1 = co, y0 = c(0, 0), y1 = c(0.4, 
            0.4) * mt, default.units = "native", gp = gpar(col = colors["highlight"], 
            lwd = 2))
    }
##############
   if (!missing(highlight3)) {
        mt = (match(highlight3$strand, c("-", "+")) - 1.5) * 2
        co = highlight3$coord
        if (is.na(mt) || !is.numeric(co))
            stop("Invalid parameter 'highlight3'.")
        grid.segments(x0 = co, x1 = co, y0 = c(0, 0), y1 = c(0.4,
            0.4) * mt, default.units = "native", gp = gpar(col = "green",
            lwd = 2))
    }
##############
   if (!missing(highlight4)) {
        mt = (match(highlight4$strand, c("-", "+")) - 1.5) * 2
        co = highlight4$coord
        if (is.na(mt) || !is.numeric(co))
            stop("Invalid parameter 'highlight4'.")
        grid.segments(x0 = co, x1 = co, y0 = c(0, 0), y1 = c(0.4,
            0.4) * mt, default.units = "native", gp = gpar(col = "green",
            lwd = 2))
    }



##############
    popViewport()
    pushViewport(viewport(layout.pos.col = 1, layout.pos.row = which(names(VP) == 
        "title")))
    grid.text(label = paste("Chr ", chr, sep = ""), x = 0.5, 
        y = 1, just = "centre", gp = gpar(cex = 1))
    if (!missing(main)) 
        grid.text(label = main, x = 0.05, y = 1, just = "centre", 
            gp = gpar(cex = 1))
    popViewport()
    if (doLegend) 
        plotAlongChromLegend(which(names(VP) == "legend"), featureColorScheme = featureColorScheme, 
            featureExclude = featureExclude)
    popViewport(2)
}

################################################################################
MyplotAlongChromLegend<-function (vpr, nr = 2, featureColorScheme = 1, featureExclude = c("chromosome", 
    "nucleotide_match", "insertion"), mainLegend, cexLegend = 0.35, 
    cexMain = 1) 
{
    endVP = FALSE
    if (missing(vpr)) {
        endVP = TRUE
        vpr = newVP(main = mainLegend, dataPanelHeight = 1, cexMain = cexMain)
    }
    formatRow = function(featColsOneRow, row) {
        strWid = convertWidth(stringWidth(rownames(featColsOneRow)), 
            "npc", valueOnly = TRUE)
        n = length(strWid)
        inbetWid = 0.2 * min(strWid)
        totWid = sum(strWid) + (n - 1) * inbetWid
        x = c(0, cumsum(strWid[-n])) + (0:(n - 1)) * inbetWid
        y = numeric(length(x))
        x = x/totWid
        strWid = strWid/totWid
        grid.rect(x = x, width = strWid, y = unit(row, "native"), 
            height = unit(1, "native") - unit(1, "mm"), just = c("left", 
                "center"), default.units = "npc", gp = do.call(gpar, 
                featColsOneRow))
        grid.text(label = rownames(featColsOneRow), x = unit(x + 
            strWid/2, "native"), y = unit(row, "native"), just = c("center", 
            "center"), gp = gpar(cex = cexLegend))
    }
    featCols = featureColors(featureColorScheme)
    featCols = featCols[!(rownames(featCols) %in% featureExclude), 
        ]
    pushViewport(viewport(layout.pos.col = 1, layout.pos.row = vpr, 
        yscale = c(0.5, nr + 0.5)))
    i = 1:nrow(featCols)
    for (r in 1:nr) formatRow(featCols[ceiling(i/nrow(featCols) * 
        nr - 1e-10) == r, ], row = nr - r + 1)
    popViewport()
    if (endVP) 
        popViewport(2)
}

################################################################################
################################################################################
