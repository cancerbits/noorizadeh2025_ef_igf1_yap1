resultsDir <- function(...) {
  paste0(config$out_root, "/", ...)
}

pToSig = function(p, ns="", threshs=c(0.05, 0.01, 0.001)) {
  sapply(p, function(pval) ifelse(pval>threshs[1], ns, ifelse(pval<=threshs[3], "***", ifelse(pval<=threshs[2], "**", "*"))) )
}

# pretty print a time difference between two proc.time() calls
time_diff <- function(start_time, end_time = NULL) {
  if (is.null(end_time)) {
    end_time <- proc.time()
  }
  dt_cpu <- lubridate::make_difftime(num = sum(end_time[c('user.self', 'sys.self')] - start_time[c('user.self', 'sys.self')]))
  dt_elapsed <- lubridate::make_difftime(num = end_time['elapsed'] - start_time['elapsed'])
  
  sprintf('Elapsed time: %1.2f %s; CPU time: %1.2f %s', 
          dt_elapsed, attr(x = dt_elapsed, which = 'units'),
          dt_cpu, attr(x = dt_cpu, which = 'units'))
}

rotateLabels <- function(angle = 45) {
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = angle, hjust=1))
}

getColors <- function(categories, pal=NULL) {	
  n <- length(categories)
  cols <- NA
  if(!is.null(pal)) {
    cols <- RColorBrewer::brewer.pal(max(n,3),pal)
  }
  else if(n == 1) cols <- "black"
  else if(n <= 8) {
    cols <- RColorBrewer::brewer.pal(max(n,3),"Set2")
    if(n < 3) cols <- cols[1:n]
  }
  else if(n <= 9) {
    cols <- RColorBrewer::brewer.pal(n,"Set1")
  }
  else if(n <= 12) {
    cols <- RColorBrewer::brewer.pal(n,"Set3")
  }
  else cols <- rainbow(n)
  return(structure(cols,names=as.character(categories)))
}

rblapply <- function(args, fun, id="id", ..., cores=1) {
  require(data.table)
  if(cores>1) {
    require(parallel)
    res <- parallel::mclapply(X=args, FUN=fun, ..., mc.cores=cores)
    names(res) <- names(args)
    res <- rbindlist(res, idcol=id, fill=T)
    res[,paste(id):=args[get(id)]] # args have been converted to numbers --> convert them back			
  } else {
    res <- rbindlist(sapply(X=args, FUN=fun, simplify=FALSE, ...), idcol=id, fill=T)
  }
  return(res)
} 
rblapplyDT <- function(dt, fun, idcol) {
  res <- apply(dt, 1, function(x) {
    res <- data.table(fun(as.list(x)))
    res[, paste(idcol):=x[idcol]]
    res
  })
  if(is.list(res)) res <- rbindlist(res)
  return(as.data.table(res))
}
repl <- function(x, a, b, by_name=F) {
  if(by_name) {
    x[names(x)==a] <- b
  } else {
    x[x==a] <- b
  }	
  x
}

msg <- function(...) {
  message(format(Sys.time(), "%Y-%m-%d %H:%M"), "\t", ...)
}

msgF <- function(...) {
  message(format(Sys.time(), "%Y-%m-%d %H:%M"), "\t", sprintf(...))
}

paste_ <- function(...) {
  paste(..., sep="_")
}

dtToGr <- function(dt, chrCol="chrom", startCol="start", endCol="end", metaCols=c()) {
  library("GenomicRanges")
  
  argList <- list()
  for(n in metaCols) {
    argList[[n]] <- dt[,get(n)]
  }
  
  argList$ranges <- IRanges(dt[,as.numeric(get(startCol))],dt[,as.numeric(get(endCol))])
  argList$seqnames <- dt[,get(chrCol)]
  
  do.call(GRanges, args=argList)
}

grToDt <- function(gr, chrCol="chrom", startCol="start", endCol="end") {		
  dt <- data.table(chrom=as.character(seqnames(gr)), start=start(gr), end=end(gr))
  setnames(dt, c(chrCol, startCol, endCol))
  dt
}

dtToDf <- function(dt, rownameCol=1) {
  df <- as.data.frame(dt)
  rownames(df) <- df[,rownameCol]
  if(is.numeric(rownameCol)) {
    df <- df[,-rownameCol,drop=FALSE] 
  }
  else {
    df <- df[,setdiff(colnames(df),rownameCol),drop=FALSE] 
  }
  df
}

abscap <- function(x, cap=0.99) {
  thresh <- quantile(abs(x), cap, na.rm=T)
  i <- !is.na(x) & abs(x)>=thresh
  x[i] <- sign(x[i]) * thresh
  x
}

absmax <- function(x) {
  x[which.max(abs(x))]
}

absmin <- function(x) {
  x[which.min(abs(x))]
}

capFirst <- function(str) paste0(toupper(substr(str, 1, 1)), substr(str, 2, nchar(str)))

loadAnnot <- function() {
  f <- list.files(file.path(config$project_root, "metadata"), pattern="samples_[^_]+.csv", full.names=T)
  dA <- rblapply(f, fread, "annot_file")
  dA[, collection:="internal"]
  
  dA[,tmp:=gsub("^S_","",sample_name)]
  dA[,tmp:=NULL]
  
  setkey(dA, "sample_name")
  dA
}

plotClusteredHeatmap <- function(d, annot, peaks, regs, samps, n, annot_colors=NULL, hclust_meth="complete", dist_meth="euclidean", focus_samples=NULL, do_cap=function(x) abscap(x, 0.95), do_scale=TRUE, fix_order_rows=FALSE, fix_order_cols=FALSE, K=2, show_rownames=F, ...) {	#cmp_short, , min_k=3 
  setkey(annot, sample_name)	
  hmData <- d[regs,samps]
  dn <- dimnames(hmData)
  if(do_scale) hmData <- t(apply(hmData, 1, scale))
  dimnames(hmData) <- dn
  hmData <- do_cap(hmData)	
  
  if(!fix_order_rows[1] | !fix_order_cols[1]) {	
    if(is.null(focus_samples)) {
      distMat <- dist(hmData, method=dist_meth)
    } else {
      distMat <- dist(hmData[,focus_samples], method=dist_meth)
    }
  }
  if(!fix_order_rows[1]) {
    rc <- fastcluster::hclust(distMat, method=hclust_meth)
  } else {
    rc <- FALSE 
    K <- NA
    hmData <- hmData[fix_order_rows,]
  }
  if(!fix_order_cols[1]) {
    cc <- fastcluster::hclust(dist(t(hmData), method=dist_meth), method=hclust_meth)
  } else {
    cc <- FALSE	
    hmData <- hmData[,fix_order_cols]
  }
  
  colAnnot <- dtToDf(annot[colnames(hmData), ]) #, FRiP, pass_qc
  rowAnnot <- dtToDf(peaks[rownames(hmData), ]) #data.frame(module_id=moduleIds)
  annotCol <- annot_colors[intersect(c(colnames(rowAnnot),colnames(colAnnot)),names(annot_colors))]
  
  p <- ComplexHeatmap::pheatmap(hmData, cellwidth=12, border_color=NA, show_rownames=show_rownames, use_raster=T, main=sprintf("%s (n = %d)", n, nrow(hmData)), cluster_col=cc, cluster_row=rc, cutree_row=K, annotation_row=rowAnnot, annotation_col=colAnnot, annotation_colors=annotCol, ...) 
  for(b in unique(colAnnot$batch)) {
    i <- colAnnot$batch==b
    p <- p + ComplexHeatmap::pheatmap(hmData[, i], cellwidth=12, border_color=NA, show_rownames=show_rownames, use_raster=T, main=paste("batch", b), cluster_col=T, cluster_row=rc, annotation_row=rowAnnot, annotation_col=colAnnot[i,,drop=F], annotation_colors=annotCol, ...) #[c(colnames(rowAnnot),colnames(colAnnot))])
  }
  if(!is.null(focus_samples)) {
    i <- intersect(colnames(hmData), focus_samples)
    p <- p + ComplexHeatmap::pheatmap(hmData[, i], cellwidth=12, border_color=NA, show_rownames=show_rownames, use_raster=T, main="selected samples", cluster_col=T, cluster_row=rc, annotation_row=rowAnnot, annotation_col=colAnnot[i,,drop=F], annotation_colors=annotCol, ...)
  }
  list(plot=p, samples=samps, clustering=rc)
}

calcTSSscore <- function(fs) {
  sapply(fs, function(f) {
    # code from PEPATAC:
    insertionsMat <- read.table(f, header=FALSE, row.names=NULL,as.is=TRUE, check.names=FALSE)
    normTSS <- insertionsMat / mean(insertionsMat[1:200,])
    TSSscore <- round(mean(normTSS[1950:2050,]),1)
    return(TSSscore)
  })
}

plotEnrichmentRes <- function(enrichment_results, top_n=3, is_sig_col="sig", clust_by=NULL, sort_plot_by="perc", rev_sort=F, term_id_col="motif_id", term_label_col="motif_label", list_col="module_id", odds_col="log2odds", pval_col="pval", padj_col="padj", freq_col="perc", freq_col_bg="perc_bg", always_include=c(), color_palettes=list(), db_col="database", n="Enrichment", only_sig=TRUE, bg_name="BG", multi_layout=TRUE, max_label_length=48) {
  
  e <- enrichment_results[ , c(term_id_col, term_label_col, is_sig_col, odds_col, pval_col, padj_col, freq_col, freq_col_bg, list_col, db_col, clust_by), with=FALSE]
  setnames(e, c("id","label","is_sig","odds","pval","padj","perc","perc_bg","list_name","db", clust_by))
  
  e[, label_short:=label]
  if(is.finite(max_label_length)) {
    e[nchar(label_short)>max_label_length, label_short := paste0(substr(label_short, 0, max_label_length-1),".."), by=label_short]
    e[nchar(label_short)<max_label_length, label_short := stringr::str_pad(label_short, max_label_length-nchar(label_short)), by=label_short]
  }
  
  #### find top-X results: ####
  
  topVars <- e	
  if(only_sig==T) topVars <- topVars[is_sig==TRUE,] 
  
  # only pick one example per motif cluster?
  if(!is.null(clust_by)) {
    topVars <- e[order(pval),][is_sig==TRUE, .SD[1,], by=c("list_name", clust_by)]
  }
  
  # rank and pick top-X:
  topVars[, rnk:=rank(pval, ties="random"), by=.(list_name)]	
  topVars <- topVars[order(list_name,rnk),][rnk<=top_n,]	
  topVars[,rn:=make.unique(id)]
  
  pData <- e[id%in%c(topVars[, id], always_include),]
  
  cols <- c("perc", "id", "label", "label_short", "odds", "padj", "list_name", "db")	
  # N.B. not actually selection the enrichement for BG, we're simply using this handle so to get only one entry per term and facet!
  bg <- unique(pData[, .(tmp=perc_bg, id, label, label_short, 0, 1, bg_name, db)])
  setnames(bg, cols)
  
  pData2 <- rbind(pData[, cols, with=F], bg)		
  pData2 <- merge(pData2, topVars[,.(term=id, sig=list_name, sig_p=pToSig(padj))], by.x="id", by.y="term", all.x=T, allow.cartesian=T)		
  pData2[is.na(sig), sig:="Others"]
  
  cols <- NULL
  if(!is.null(n) & paste0("",n)!="") cols <- color_palettes[[n]]
  cols <- c(cols, getColors(pData2[,unique(list_name)]))
  cols[bg_name] <- "grey"
  
  ub <- pData2[,max(perc*100,na.rm=T)]
  
  pData2[, sort_col:=get(sort_plot_by)]
  if(rev_sort==T) pData2[, sort_col:=-sort_col]
  
  p <- ggplot(pData2, aes(y=perc*100, x=reorder(sprintf("%s", label_short), sort_col))) + geom_bar(stat="identity", width=0.75, position=position_dodge(width=0.9), aes(fill=list_name)) + theme(legend.position="top",legend.title=element_blank()) + coord_flip() + theme(panel.border=element_rect(colour="black", fill=NA)) + xlab(NULL) + ylab("Overlap (%)") + scale_fill_manual(values=cols) + scale_y_continuous(expand=c(0,0), limits=c(0, ub*1.2))
  p <- p + geom_text(aes(y=ub*1.19, label=sig_p), hjust=1, data=pData2[list_name==sig,])	

  if(multi_layout) {
    p <- p + facet_wrap(sprintf("Enriched in: %s",sig)~db, ncol=pData2[,length(unique(db))], scales="free", drop=F) 
  } else {
    if(pData2[,length(unique(db))]>1) {
      p <- p + ggforce::facet_col(db~sprintf("Enriched in: %s",sig), scales="free_y", space="free", drop=F) 
    } else {
      p <- p + ggforce::facet_col(.~sprintf("Enriched in: %s",sig), scales="free_y", space="free", drop=F) + ggtitle(pData2[,unique(db)])
    }
  }
  
  return(list(plot=p, results=e, top_vars=topVars))
}

# modified function from fgsea:
plotEnrichmentMod <- function (pathway, stats, mn, rn, padj, nes, dar_dir, gseaParam = 1, ticksSize = 0.2, label_gene = NULL)
{
  pw <- pathway
  rnk <- rank(-stats)
  ord <- order(rnk)
  statsAdj <- stats[ord]
  statsAdj <- sign(statsAdj) * (abs(statsAdj)^gseaParam)
  statsAdj <- statsAdj/max(abs(statsAdj))
  pathway <- unname(as.vector(na.omit(match(pathway, names(statsAdj)))))
  pathway <- sort(pathway)
  gseaRes <- fgsea::calcGseaStat(statsAdj, selectedStats = pathway,
                                 returnAllExtremes = TRUE)
  bottoms <- gseaRes$bottoms
  tops <- gseaRes$tops
  n <- length(statsAdj)
  xs <- as.vector(rbind(pathway - 1, pathway))
  ys <- as.vector(rbind(bottoms, tops))
  toPlot <- data.frame(x = c(0, xs, n + 1), y = c(0, ys, 0))
  diff <- (max(tops) - min(bottoms))/8
  x = y = NULL
  g1 <- ggplot(toPlot, aes(x = x, y = y)) + 
        ylab("Enrichment score") + xlab(NULL) +
        ggtitle(sprintf("%s: %s, padj=%.3f, NES=%.1f", rn, mn, padj, nes)) +
        geom_point(color = "green",size = 0.1) + 
        geom_hline(yintercept = max(tops), colour = "red", linetype = "dashed") + 
        geom_hline(yintercept = min(bottoms), colour = "red", linetype = "dashed") + 
        geom_hline(yintercept = 0, colour = "black") + 
        geom_line(color = "green") +
        geom_segment(data = data.frame(x = pathway), mapping = aes(x = x, y = -diff/2, xend = x, yend = diff/2), linewidth = ticksSize) +
        theme(panel.border = element_blank(), panel.grid.minor = element_blank(), axis.text.x=element_blank()) 
  
  dx <- data.table(x=1:length(statsAdj), v=-statsAdj, sig_dir=as.factor(dar_dir[names(statsAdj)]), in_pathway=names(statsAdj)%in%pw)
  g2 <- ggplot(dx, aes(x = x, y = v, color=sig_dir)) + 
        geom_line(color="black") +
        geom_hline(yintercept=0, linetype="dashed") + 
        geom_point(data=dx[sig_dir!=0,], shape="+", size=1.5) +
        xlab("Peaks ranked by DESeq2 test statistic") + ylab("DEseq2 stat") + scale_color_manual(values=colorPalettes$sig_dir, guide=F) 
  
  if(!is.null(label_gene)) {
    i <- which(names(statsAdj)==toupper(label_gene))
    g2 <- g2 + annotate(geom="text", x=i, y=statsAdj[i], label=label_gene)
  }
  
  patchwork::wrap_plots(g1, g2, ncol=1, tag_level="new")
}

makeDirNames <- function(sig_dir, cmp) {
  cmp <- gsub("_sig_dir", "", cmp)
  grpNames <- strsplit(cmp, "_vs_")[[1]]
  dirNames <- c(
    sprintf("%s < %s, n = %d", grpNames[1], grpNames[2], sum(sig_dir=="-1")),
    sprintf("%s ~ %s, n = %d", grpNames[1], grpNames[2], sum(sig_dir=="0")),
    sprintf("%s > %s, n = %d", grpNames[1], grpNames[2], sum(sig_dir=="1"))
  )
  dirNames2 <- c(
    sprintf("%s < %s", grpNames[1], grpNames[2]),
    sprintf("%s ~ %s", grpNames[1], grpNames[2]),
    sprintf("%s > %s", grpNames[1], grpNames[2])
  )
  return(list(
    dir_names = dirNames,
    dir_names_nonum = dirNames2,
    grp_names = grpNames,
    replaced = dirNames[as.numeric(sig_dir)+2],
    replaced_nonum = dirNames2[as.numeric(sig_dir)+2]
  ))
}
putInBin <- function(numbers, bins) {
  sapply(numbers, function(x) {
    which(sapply(bins, function(b) x>=b[1] & x<=b[2]))[1]
  })
}