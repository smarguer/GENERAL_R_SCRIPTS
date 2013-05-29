
plotGENE=function (gff, chr, xlim, strand, vpr, featureColorScheme = 1, 
    featureExclude = c("chromosome", "nucleotide_match", "insertion"), 
    featureNoLabel = c("uORF", "CDS"), ...) 
{
    pushViewport(dataViewport(xData = xlim, yscale = c(-1.2, 
        1.2), extension = 0, clip = "on", layout.pos.col = 1, 
        layout.pos.row = vpr))
    stopifnot(all(gff[, "start"] <= gff[, "end"]))
    sel = which(gff[, "chr"] == chr & gff[, "strand"] == strand & 
        gff[, "start"] <= xlim[2] & gff[, "end"] >= xlim[1])
    stopifnot(length(strand) == 1, strand %in% c("+", "-"))
    geneName = gff[sel, "gene"]
    featName = gff[sel, "Name"]
    featName[!is.na(geneName)] = geneName[!is.na(geneName)]
    feature = as.character(gff[sel, "feature"])
    featsp = split(seq(along = sel), feature)
    whnames = integer(0)
    featsp = featsp[!(names(featsp) %in% featureExclude)]
    wh = ("gene" == names(featsp))
    if (any(wh)) {
        i = featsp[["gene"]]
        s = sel[i]
        grid.segments(x0 = gff$start[s], x1 = gff$end[s], y0 = 0, 
            y1 = 0, default.units = "native", gp = gpar(col = "#a0a0a0"))
        whnames = i
        featsp = featsp[!wh]
    }
    wh = ("intron" == names(featsp))
    if (any(wh)) {
        i = featsp[["intron"]]
        s = sel[i]
        mid = (gff$start[s] + gff$end[s])/2
        wid = (gff$end[s] - gff$start[s])/2
        for (z in c(-1, 1)) grid.segments(x0 = mid, x1 = mid + 
            z * wid, y0 = 1.2 * c(`+` = 1, `-` = -1)[strand], 
            y1 = 0.95 * c(`+` = 1, `-` = -1)[strand], default.units = "native", 
            gp = gpar(col = "black"))
        featsp = featsp[!wh]
    }
    featCols = featureColors(featureColorScheme)
    whm = names(featsp) %in% rownames(featCols)
    if (!all(whm)) 
        warning("Don't know how to handle feature of type(s) '", 
            paste(names(featsp)[!whm], collapse = ", "), "' in gff.", 
            sep = "")
    sfeatsp = featsp[rownames(featCols)]
    ll = listLen(sfeatsp)
    if (any(ll > 0)) {
        i = unlist(sfeatsp)
        gp = gpar(col = rep(featCols$col, ll), fill = rep(featCols$fill, 
            ll))
        s = sel[i]
        grid.rect(x = gff$start[s], y = 0, width = gff$end[s] - 
            gff$start[s], height = 2, default.units = "native", 
            just = c("left", "center"), gp = gp)
        whnames = c(whnames, unlist(sfeatsp[!(names(sfeatsp) %in% 
            featureNoLabel)]))
    }
    if (!all(tolower(featureNoLabel) == "all") && (length(whnames) > 
        0)) {
        bindingRegexpr = "binding.?site.*$"
        isBindingSite = (regexpr(bindingRegexpr, featName[whnames]) > 
            0)
        if (any(isBindingSite)) {
            featName[whnames] = gsub(bindingRegexpr, "bs", featName[whnames])
        }
        whnames = whnames[isBindingSite | !duplicated(featName[whnames])]
        txtcex = 0.6
        txtdy = 0.7
        s = sel[whnames]
        txtx = (gff$start[s] + gff$end[s])/2
        txty = numeric(length(s))
        ord = order(txtx)
        whnames = whnames[ord]
        s = s[ord]
        txtx = txtx[ord]
        strw = convertWidth(stringWidth(featName[whnames]), "native", 
            valueOnly = TRUE) * txtcex
        rightB = txtx[1] + 0.5 * strw[1]
        doText = rep(TRUE, length(whnames))
        if (length(whnames) > 1) {
            for (k in 2:length(whnames)) {
                leftB = txtx[k] - 0.5 * strw[k]
                if (leftB > rightB) {
                  rightB = txtx[k] + 0.5 * strw[k]
                }
                else {
                  if (!any(txty[k - (1:2)] == txtdy)) {
                    txty[k] = txtdy
                  }
                  else {
                    if (!any(txty[k - (1:2)] == -txtdy)) {
                      txty[k] = -txtdy
                    }
                    else {
                      doText[k] = FALSE
                    }
                  }
                }
            }
        }
        grid.text(label = featName[whnames][doText], x = txtx[doText], 
            y = txty[doText], gp = gpar(cex = txtcex), default.units = "native")
    }
    popViewport()
}
