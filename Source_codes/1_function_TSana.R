## calculate each isoform relative ration in gene expression
run_ratioIso <- function(iso_geneInfor,rawexpmat,type=c("ratio","geneSUM")){
  colnames(iso_geneInfor) <- c("iso","geneID")
  rownames(iso_geneInfor) <- iso_geneInfor$iso
  geneSummat <- rowsum(rawexpmat,group =iso_geneInfor$geneID,reorder = F )
  if(type=="ratio"){
    geneSummat <- geneSummat[iso_geneInfor$geneID,]
    if(is.null(nrow(rawexpmat))){
      rawexpmat <- rawexpmat[iso_geneInfor$iso]
    }else{
      rawexpmat <- rawexpmat[iso_geneInfor$iso,]
    }
    
    ratioMat <- rawexpmat/geneSummat
    ratioMat[is.na(ratioMat)] <- 0
    return(ratioMat)
  }else if(type=="geneSUM"){
    return(geneSummat)
  }
  
}
## mean replicate expression
isoMat.meanRep <- function(isoRatio_mat,group_rep){
  Meanrep_isomat <- t(isoRatio_mat)
  Meanrep_isomat <- rowsum(Meanrep_isomat,group = group_rep,reorder = F)
  Meanrep_isomat <- t(Meanrep_isomat)
  rep_num <- table(group_rep)[colnames(Meanrep_isomat)]
  rep_num <- matrix(rep_num,nrow = nrow(Meanrep_isomat),ncol = length(rep_num),byrow = T)
  Meanrep_isomat <- Meanrep_isomat/rep_num;rm(rep_num)
  return(Meanrep_isomat)
}
## extract isoform id by ration threhold

## Zscores
zscore_cal <- function(x){(x-mean(x))/sd(x)}

# plot expression ---------------------------------------------------------
bed_to_grange <- function(bedMat,strand_if = F){
  bedMat <- data.frame(bedMat,stringsAsFactors = F)
  bedMat[,2] <- as.numeric(bedMat[,2])
  bedMat[,3] <- as.numeric(bedMat[,3])
  if(strand_if){
    tmp_gr <- GRanges(seqnames =  Rle(bedMat[,1]), ranges = IRanges(start = bedMat[,2],width =bedMat[,3] - bedMat[,2]+1),
                      strand = Rle(bedMat[,4]))
    bedMat <- bedMat[,-4]
  }else{
    tmp_gr <- GRanges(seqnames =  Rle(bedMat[,1]), ranges = IRanges(start = bedMat[,2],width =bedMat[,3] - bedMat[,2]+1),
                      strand = Rle(strand(rep("*",nrow(bedMat)))))
  }
  
  if(ncol(bedMat) > 3){
    mcols(tmp_gr) <- bedMat[,4:ncol(bedMat)]
    colnames(mcols(tmp_gr)) <- colnames(bedMat)[4:ncol(bedMat)]
  }
  
  tmp_gr
}
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  require(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- plyr::ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}



plot_expFile <- function(iso_plotID,time_set,iso_genepair,plot_expMat,ifplot_geneExp=T,
                         ifall_iso=T,plot_ratio=F,ylab="TPM",ifoutput_mat=F,if_replot=F,errorbar=T,crossMat=NULL,mainTitle=NULL,color_line=NULL,
                         remove_iso=NULL){
  if(length(time_set) != ncol(plot_expMat)){
    return(NULL)
  }
  iso_plotID <- as.character(iso_plotID)
  colnames(iso_genepair) <- c("iso","geneID")
  rownames(iso_genepair) <- iso_genepair$iso
  
  gene_name <- iso_genepair[iso_plotID[1],"geneID"]
  allisoID <- iso_genepair[which(iso_genepair$geneID ==gene_name) ,"iso"]
  oneGene_exp <- apply(plot_expMat[allisoID,],2,sum)
  
  if(ifall_iso){
    iso_allId <- allisoID
  }else{
    iso_allId <- iso_plotID
  }
  if(length(remove_iso) >0){
    iso_allId <- setdiff(iso_allId,remove_iso)
  }
  if(is.null(color_line)){
    color_line <- rainbow(length(iso_allId))
  }
 

  
  names(color_line) <-iso_allId
  line_type <- c(rep(1,length(iso_plotID)),rep(2,length(setdiff(iso_allId,iso_plotID)) ))
  names(line_type) <- c(iso_plotID,setdiff(iso_allId,iso_plotID))
  ###
  if(ifplot_geneExp){
    color_line <- c( "black",color_line)
    gene_name <- paste0("gene: ",gene_name)
    names(color_line)[1] <- gene_name
    line_type <- c(1,line_type)
    names(line_type)[1] <- gene_name
    one_expmat <- rbind(oneGene_exp,
                        plot_expMat[iso_allId,])
    rownames(one_expmat)[1] <- gene_name
  }else{
    one_expmat <- plot_expMat[iso_allId,]  }
  

  
  if(plot_ratio){
    one_expmat <- rbind(oneGene_exp,one_expmat)
    one_expmat <- apply(one_expmat,2,function(cc){cc[-1]/cc[1]})
    one_expmat[is.na(one_expmat)] <- 0
    # crossMat <- NULL
  }
  
  ### plot
  one_expmat <- data.frame(one_expmat,stringsAsFactors = F)
  one_expmat$id <- rownames(one_expmat)
  time_setNew <- rep(time_set,each=nrow(one_expmat))
  one_expmat <- reshape2::melt(one_expmat)
  one_expmat$variable <- time_setNew
  colnames(one_expmat) <- c("id","times","expValue")
  allrep_mat <- one_expmat
  try(one_expmat <- summarySE(one_expmat, measurevar="expValue", groupvars=c("id","times")))
  if(ifoutput_mat){
    return(list(one_expmat=one_expmat,color_line=color_line,line_type=line_type))
  }
  if(if_replot){
    allrep_matplot <- allrep_mat
  }else{
    allrep_matplot <- one_expmat
  }
  # ggplot(one_expmat, aes(x=times, y=expValue,colour=id,group=id)) + 
  #   geom_errorbar(aes(ymin=expValue-se, ymax=expValue+se), width=.1) +
  #   geom_line(aes(linetype=id)) +
  #   geom_point(data =allrep_mat ,aes(x=times, y=expValue,colour=id,group=id))+
  #   scale_color_manual(values=color_line)+
  #   scale_linetype_manual(values =line_type)+
  #   ylab(label =ylab)
  # 
  
  if(errorbar){
    try(plotobj <- ggplot(one_expmat, aes(x=times, y=expValue,colour=id,group=id)) + 
          geom_errorbar(aes(ymin=expValue-se, ymax=expValue+se), width=.1) +
          geom_line(aes(linetype=id)) +
          geom_point(data =allrep_matplot ,aes(x=times, y=expValue,colour=id,group=id)) +
          scale_color_manual(values=color_line)+
          scale_linetype_manual(values =line_type)+
          ylab(label =ylab) )
  }else{
    try(plotobj <- ggplot(one_expmat, aes(x=times, y=expValue,colour=id,group=id)) + 
          geom_line(aes(linetype=id)) +
          geom_point(data =allrep_matplot ,aes(x=times, y=expValue,colour=id,group=id)) +
          scale_color_manual(values=color_line)+
          scale_linetype_manual(values =line_type)+
          ylab(label =ylab))
  }
  if(!is.null(crossMat)){
    plotobj <- plotobj + 
      annotate("pointrange", x =crossMat$x.value,
               y = crossMat$y.value, ymin = 0, ymax = 0,
               colour = "black", size = 0.5)
  }
  plotobj <-plotobj+labs(title = mainTitle)
  return(plotobj)
}


plot_expFilehtml <- function(iso_plotID,time_set,iso_genepair,plot_expMat,
                             ylab="TPM",xlab="Time",errorbar=T,crossMat=NULL,mainTitle=NULL){
  if(length(time_set) != ncol(plot_expMat)){
    return(NULL)
  }
  iso_plotID <- sort(as.character(iso_plotID))
  colnames(iso_genepair) <- c("iso","geneID")
  rownames(iso_genepair) <- iso_genepair$iso
  
  gene_name <- iso_genepair[iso_plotID[1],"geneID"]
  allisoID <- iso_genepair[iso_genepair$geneID ==gene_name ,"iso"]
  oneGene_exp <- apply(plot_expMat[allisoID,],2,sum)
  iso_allId <- iso_plotID
  color_line <- c( "#1660a7" ,"#cd0c18") #rainbow(length(iso_allId))
  names(color_line) <-iso_allId
  line_type <- c(rep(1,length(iso_plotID)),rep(2,length(setdiff(iso_allId,iso_plotID)) ))
  names(line_type) <- c(iso_plotID,setdiff(iso_allId,iso_plotID))
  ###
  one_expmat <- plot_expMat[iso_allId,]
  ### plot
  one_expmat <- data.frame(one_expmat,stringsAsFactors = F)
  one_expmat$id <- rownames(one_expmat)
  time_setNew <- rep(time_set,each=nrow(one_expmat))
  one_expmat <- reshape2::melt(one_expmat)
  one_expmat$variable <- time_setNew
  colnames(one_expmat) <- c("id","times","expValue")
  try(one_expmat <- summarySE(one_expmat, measurevar="expValue", groupvars=c("id","times")))
  
  if(errorbar){
    try(plotobj <- plot_ly(data = one_expmat[one_expmat$id == names(color_line)[1],], x = ~times, y = ~expValue, type = 'scatter', mode = 'lines+markers',
                           name = names(color_line)[1],color=I(color_line[1]), line = list( width = 4),marker=list(size=10),
                           error_y = ~list(array = se)) %>%
          add_trace(data =  one_expmat[one_expmat$id == names(color_line)[2],], 
                    name = names(color_line)[2], color=I(color_line[2])) %>%
          layout( title =mainTitle,
                  paper_bgcolor='rgb(255,255,255)', plot_bgcolor='#EBEBEB',
                  xaxis = list(title=xlab,
                               gridcolor = 'rgb(255,255,255)',
                               showgrid = TRUE,
                               showline = FALSE,
                               showticklabels = TRUE,
                               tickcolor = 'rgb(127,127,127)',
                               ticks = 'outside',
                               zeroline = FALSE),
                  yaxis = list(title = ylab,
                               gridcolor = 'rgb(255,255,255)',
                               showgrid = TRUE,
                               showline = FALSE,
                               showticklabels = TRUE,
                               tickcolor = 'rgb(127,127,127)',
                               ticks = 'outside',
                               zeroline = FALSE)
          ))
  }else{
    try(plotobj <- plot_ly(data = one_expmat[one_expmat$id == names(color_line)[1],], x = ~times, y = ~expValue, type = 'scatter', mode = 'lines+markers',
                           name = names(color_line)[1],color=I(color_line[1]),line = list( width = 4),marker=list(size=10)  ) %>%
          add_trace(data =  one_expmat[one_expmat$id == names(color_line)[2],], 
                    name = names(color_line)[2], color=I(color_line[2])) %>%
          layout( title =mainTitle,
                  paper_bgcolor='rgb(255,255,255)', plot_bgcolor='#EBEBEB',
                  xaxis = list(gridcolor = 'rgb(255,255,255)',
                               showgrid = TRUE,
                               showline = FALSE,
                               showticklabels = TRUE,
                               tickcolor = 'rgb(127,127,127)',
                               ticks = 'outside',
                               zeroline = FALSE),
                  yaxis = list(title = ylab,
                               gridcolor = 'rgb(255,255,255)',
                               showgrid = TRUE,
                               showline = FALSE,
                               showticklabels = TRUE,
                               tickcolor = 'rgb(127,127,127)',
                               ticks = 'outside',
                               zeroline = FALSE)
          ))
  }
  if(!is.null(crossMat)){
    plotobj <- plotobj %>%  add_trace(data = crossMat,x = ~x.value,y = ~y.value, name = 'cross point', 
                                      mode = 'markers',error_y=NULL,marker=list(size = 12,color = "black"),line = list( width = 0))  }
  return(plotobj)
}



plot_pairwiseTS <- function(iso_plotID,condition_set,iso_genepair,plot_expMat,ylab="TPM"){
  if(length(condition_set) != ncol(plot_expMat)){
    return(NULL)
  }
  iso_plotID <- as.character(iso_plotID)
  colnames(iso_genepair) <- c("iso","geneID")
  rownames(iso_genepair) <- iso_genepair$iso
  
  gene_name <- iso_genepair[iso_plotID[1],"geneID"]
  allisoID <- iso_genepair[which(iso_genepair$geneID ==gene_name) ,"iso"]
  oneGene_exp <- apply(plot_expMat[allisoID,],2,sum)
  
  one_expmat <- plot_expMat[iso_plotID,]
  one_expmat <- data.frame(one_expmat,stringsAsFactors = F)
  one_expmat <- isoMat.meanRep(one_expmat,condition_set)
  one_expmat <- data.table::melt(one_expmat)
  colnames(one_expmat) <- c("id","condition","expValue")
  ggplot(one_expmat,aes(x=condition,y=expValue,fill=id)) + 
    geom_bar(stat = "identity",position = "dodge") +
    ylab(label =ylab)
 
 
}
#extrat information form sqlite
extract_sqlite <- function(sqliteName,isoID=NULL,col_id=NULL,if_exp=F,extract_wholetable=NULL){
  sqliteNb <- dbConnect(SQLite(), sqliteName)
  if(!is.null(extract_wholetable)){
    # idname_mat <- dbGetQuery(sqliteNb, "select distinct isoID, geneID from iso_Structure")
    idname_mat <- dbReadTable(sqliteNb,extract_wholetable)
    return(idname_mat)
  }
  if(if_exp){
    if(is.null(col_id)){
      col_id <- "*"
    }else{
      col_id <- do.call("paste",c(as.list(c("isoID",col_id)),sep=", "))
    }
    if(is.null(isoID)){
      query_str <- paste0("select ",col_id," from iso_abundance")
    }else{
      isoID <- lapply(isoID,function(ii){paste0("'",ii,"'")})
      isoID <- do.call("paste",c(isoID,sep=", "))
      query_str <- paste0("select ",col_id," from iso_abundance where isoID in (",
                          isoID,")")
      
    }
    oneExpmat  <-   dbGetQuery(sqliteNb, query_str)
    rownames(oneExpmat) <- oneExpmat$isoID
    dbDisconnect(sqliteNb)
    return(oneExpmat[,-1])
  }
  query_str <- paste0("select seqnames, start, end, width, strand, type, isoID from iso_Structure where isoID in ('",
                      isoID,"')")
  isoform_infor <-  dbGetQuery(sqliteNb, query_str)
  dbDisconnect(sqliteNb)
  return(isoform_infor)
  
}
paste_vector <- function(string_vector,sep=";"){
  string_vector <-  base::Reduce(paste0,paste0(string_vector,sep),)
  sub(paste0(sep,"$"),"",string_vector)
} 
extract_SNPsqlite <- function(sqliteName,tableName=NULL,fields_id=NULL,query_field=NULL,
                              match_ID=NULL,match_type="=",extract_wholetable=FALSE){
  sqliteNb <- dbConnect(SQLite(), sqliteName)
  if(is.null(tableName)){
    dbDisconnect(sqliteNb)
    return(tableName)
  }
  if(extract_wholetable){
    # idname_mat <- dbGetQuery(sqliteNb, "select distinct isoID, geneID from iso_Structure")
    tablemat <- dbReadTable(sqliteNb,tableName)
    dbDisconnect(sqliteNb)
    return(tablemat)
    
  }
  if(is.null(fields_id)| fields_id[1] == "*"){
    fields_id <- "*"
  }else{
    fields_id <- do.call("paste",c(as.list(fields_id),sep=", "))
  }
  if(is.null(query_field) | is.null(match_ID)){
    query_str <- paste0("select ",fields_id," from ",tableName)
    one_extractMat  <-   dbGetQuery(sqliteNb, query_str)
  }else{
    if(length(match_ID) < 100){
      query_str <- paste0("select ",fields_id," from ",tableName," where ",query_field," ",match_type," :ID")
      one_extractMat  <-   dbGetQuery(sqliteNb, query_str, params=list(ID= match_ID))
    }else{
      match_ID <- paste0("'",match_ID,"'")
      match_ID <- paste_vector(match_ID,sep = ",")
      query_str <- paste0("select ",fields_id," from ",tableName," where ",query_field," in (",match_ID,")")
      one_extractMat  <-   dbGetQuery(sqliteNb, query_str)
    }
    
  }
  dbDisconnect(sqliteNb)
  return(one_extractMat)
}


# extract isoform sturcture
isoform_plotstructure_old <- function(structureMat,type_col=NULL){
  if(is.null(type_col)){
    type_col <- c("blue","red","purple")
    names(type_col) <- c("5' utr","coding_region","3' utr" )
  }
  xlim_x <- range(structureMat[,c("start","end")])
  if(structureMat$strand[1] == "-"){
    structureMat[,2:3] <-as.matrix( structureMat[,3:2])
  }
  structureMat <- structureMat[order(structureMat$isoID,structureMat$start),]
  structureMat$yend <- rep(1:length(unique(structureMat$isoID)),table(structureMat$isoID))
  intro_mat <- structureMat[structureMat$type == "intron",]
  other_mat <- structureMat[structureMat$type != "intron",]
  fill_col <- type_col[other_mat$type]
  
  x_scaleseq <- seq(xlim_x[1],xlim_x[2],(xlim_x[2]-xlim_x[1])/2) 
  x_scaleseq_lab <- genomePos_scalename(x_scaleseq)
  
  ggplot(other_mat, aes(xmin = start, xmax = end, y = isoID, fill = type,colour=type)) +
    geom_gene_arrow(arrowhead_height = unit(5, "mm"), arrowhead_width = unit(1.1, "mm"),
                    arrow_body_height = unit(5, "mm")) +  # geom_gene_label(align = "left")
    # geom_segment(data =intro_mat,aes(x=start,xend=end,y=isoID,yend=isoID) )
    annotate("segment",x=intro_mat$start,xend = intro_mat$end,
             y=intro_mat$yend,yend =intro_mat$yend,size=1.5,colour="grey")+
    annotate("text", x =  mean(xlim_x),y = unique(intro_mat$yend),hjust =0.5,vjust =-2.9,angle = 0,
             label = unique(intro_mat$isoID),
             size = 4,colour  = "black") +
    scale_fill_manual(values = fill_col)+
    scale_color_manual(values=fill_col )+
    xlim(xlim_x)+
    scale_x_continuous(breaks = x_scaleseq,labels = x_scaleseq_lab)+
    theme( panel.background=element_blank(),axis.text.y  = element_blank(),axis.title.x=element_blank(),
           axis.title.y=element_blank(),axis.ticks.y = element_blank()) 
  
}

# structureMat <- plot_iso_structure
isoform_plotstructure <- function(structureMat){
  xlim_x <- range(structureMat[,c("start","end")])
  # require(ggbio)
  structureMat <- bed_to_grange(structureMat[,-4],strand_if = F)
  strand(structureMat) <- structureMat$strand
  structureMat <- structureMat[,-1]
  intro_obj <- structureMat[structureMat$type == "intron"]
  structureMat <- structureMat[structureMat$type != "intron"]
  structureMat <- split(structureMat,structureMat$isoID)
  intro_obj <- split(intro_obj,intro_obj$isoID)
  # x_scaleseq <- seq(xlim_x[1],xlim_x[2],(xlim_x[2]-xlim_x[1])/2) 
  # x_scaleseq_lab <- genomePos_scalename(x_scaleseq)
  ggplot()+ ggbio::geom_segment(structureMat,stat = "identity",aes(color=type,x=start,y=isoID),size=5) + 
    ggbio::geom_arrow(intro_obj,stat = "identity",aes(x=start,y=isoID),arrow.rate =0.02,angle=45) +
    ggplot2::annotate("text", x =  mean(xlim_x),y = names(structureMat),hjust =0.5,vjust =-2.9,angle = 0,
                      label = names(structureMat),
                      size = 4,colour  = "black")   + scale_color_brewer(palette = "Set1") +
    theme( panel.background=element_blank(),axis.text.y  = element_blank(),axis.title.x=element_blank(),
           axis.title.y=element_blank(),axis.ticks.y = element_blank())
  
  # ggplot()+   ggbio::geom_chevron(intro_obj,stat = "stepping",aes(group=isoID),offset=0.05)  +
  #   ggbio::geom_segment(structureMat,stat = "stepping",aes(color=type,group=isoID),size=5) +  
  #   ggplot2::annotate("text", x =  mean(xlim_x),y = 1:length(structureMat), hjust =0.5,vjust =-3,angle = 0,
  #                     label = names(structureMat),
  #                     size = 4,colour  = "black")   + 
  #   theme( panel.background=element_blank(),axis.text.y  = element_blank(),axis.title.x=element_blank(),
  #                                                         axis.title.y=element_blank(),axis.ticks.y = element_blank()) 
    # scale_x_continuous(breaks = x_scaleseq,labels = x_scaleseq_lab)

  # 
  # ggplot()+ ggbio::geom_arrowrect(structureMat,stat = "stepping",aes(fill=type,y=start,group=isoID)) + 
  #   ggbio::geom_chevron(intro_obj,stat = "stepping",aes(group=isoID))
  
}


genomePos_scalename <- function(aa){
  max_value <- max(aa)
  if(max_value < 1e4){
    paste0(aa,"bp")
  }else if(max_value < 1e7 ){
    paste0(round(aa/1e3,3),"kb")
  }else{
    paste0(round(aa/1e6,3),"Mb")
  }
}


# extract isoform protein sturcture and plot
isoform_proteinplot <- function(isoform_ids,all_isopepMat,strand="+"){
  domain_mat <-  all_isopepMat[all_isopepMat[,1] %in% isoform_ids,c(1,3:4,6,5,2)]
  if(nrow(domain_mat) == 0){
    return(ggplot() + theme_void())
  }
  y_sort <- 1:length(isoform_ids)
  names(y_sort) <-  sort(isoform_ids)
  colnames(domain_mat) <- c("isoID","start","end","domin","pro_len","pfam")
  protein_mat <- data.frame(isoID=domain_mat$isoID,start=1,end=domain_mat$pro_len,stringsAsFactors = F)
  protein_mat <- unique(protein_mat)
  xlim_x <- max(domain_mat$pro_len,na.rm = T)
  ## if strand == "-"
  if(strand=="-"){
    domain_mat[,2:3] <- xlim_x- as.matrix(domain_mat[,3:2])+1
    protein_mat[,2:3] <- xlim_x- as.matrix(protein_mat[,3:2])+1
  }
  xlim_x <- c(1,xlim_x)
  # xlim_x <- range(protein_mat[,c("start","end")])
  
  #rep(y_sort[unique(domain_mat$isoID)],table(domain_mat$isoID))
  ## if no domain by predicted
  domain_mat <- domain_mat[!is.na(domain_mat$start),]
  if(nrow(domain_mat) == 0 & nrow(protein_mat) != 0){
    plotobj_pro <- ggplot(protein_mat)+
      geom_segment(aes(x = start, xend = end, y = isoID,yend=isoID,size=0.3))+
      theme( panel.background=element_blank(),axis.text.y  = element_blank(),axis.title.x=element_blank(),
             axis.title.y=element_blank(),axis.ticks.y = element_blank(),axis.ticks =element_blank(),
             axis.text.x = element_blank(),legend.position = "none" )
    return(plotobj_pro)
    
  }
  domain_mat$yend <- y_sort[domain_mat$isoID] 
  unique_domain <- unique(domain_mat$domin)
  if(length(unique_domain) < 3){
    domain_col <- rainbow(length(unique_domain))
  }else{
    domain_col <- RColorBrewer::brewer.pal(length(unique_domain),"Set1")
  }
  names(domain_col) <- unique_domain
  domain_mat$col <- domain_col[domain_mat$domin]
  protein_mat <- bed_to_grange(protein_mat)
  tt_doma <- bed_to_grange(domain_mat)
  nodomain_mat  <- GenomicRanges::setdiff(protein_mat,tt_doma)
  nodomain_mat <- as.data.frame(nodomain_mat,stringsAsFactors=F)[,1:3]
  domain_mat$width <- domain_mat$end -domain_mat$start +1
  domain_mat <- domain_mat[order(domain_mat$isoID,domain_mat$start),]
  # if only one trancript codeing protein
  add_id <- setdiff(isoform_ids,nodomain_mat$seqnames)
  if(length(add_id) != 0){
    nodomain_mat <- rbind(nodomain_mat,data.frame(seqnames=add_id,start=1,end=1,stringsAsFactors = F))
  }
  # GenomicRanges::gaps(tt_doma,start =protein_mat$start ,end =protein_mat$end)
  domain_mat$isoID <- factor(domain_mat$isoID,levels = sort(isoform_ids))
  plotobj_pro <- ggplot(domain_mat)+
    # geom_segment(aes(x = start, xend = end, y = isoID,yend=isoID,color= domain_mat$domin,size=2))
    geom_segment(aes(x = start, xend = end, y = isoID,yend=isoID, colour=domin,size=2)) +
    geom_segment(data =nodomain_mat ,aes(x = start, xend = end, y = seqnames,yend=seqnames,size=1)) +
    annotate("text", x =  domain_mat$start + domain_mat$width*0.1 ,y = domain_mat$yend,hjust =0.2,vjust =-1.3,angle = 0,
             label = domain_mat$domin,
             size = 3,colour  = "black") +
    annotate("text", x =  mean(xlim_x),y = isoform_ids,hjust =0.5,vjust =-2.9,angle = 0,
             label = paste0(isoform_ids," Protein"),
             size = 4,colour  = "black")  +
    # theme( panel.background=element_blank(),axis.text.y  = element_blank(),axis.title.x=element_blank(),
    #        axis.title.y=element_blank(),axis.ticks.y = element_blank()) + guides(size=FALSE) 
    # scale_fill_manual(values = domain_col)+scale_color_manual(values=fill_col )
    theme( panel.background=element_blank(),axis.text.y  = element_blank(),axis.title.x=element_blank(),
           axis.title.y=element_blank(),axis.ticks.y = element_blank(),axis.ticks =element_blank(),
           axis.text.x = element_blank(),legend.position = "none" )
  plotobj_pro
  
}

# venn plot ---------------------------------------------------------------

extract_TSinfor <- function(onemat,type="isopair"){
  if(type == "isopair"){
    iso1_id <-pmin(as.character(onemat$iso1),as.character(onemat$iso2))
    iso2_id <-pmax(as.character(onemat$iso1),as.character(onemat$iso2))
    return(unique(paste0(iso1_id,"&",iso2_id)))
  }
  if(type == "gene"){
    gene_id <- onemat$geneID #isoform_geneid[c(as.character(onemat$iso1),as.character(onemat$iso2)),1]
    return(unique(gene_id))
  }
  
}


venndiagramNew <- function (x, filename, height = 3000, width = 3000, resolution = 500, 
                            imagetype = "tiff", units = "px", compression = "lzw", na = "stop", 
                            main = NULL, sub = NULL, main.pos = c(0.5, 1.05), main.fontface = "plain", 
                            main.fontfamily = "serif", main.col = "black", main.cex = 1, 
                            main.just = c(0.5, 1), sub.pos = c(0.5, 1.05), sub.fontface = "plain", 
                            sub.fontfamily = "serif", sub.col = "black", sub.cex = 1, 
                            sub.just = c(0.5, 1), category.names = names(x), force.unique = TRUE, 
                            print.mode = "raw", sigdigs = 3, direct.area = FALSE, area.vector = 0, 
                            hyper.test = FALSE, total.population = NULL, lower.tail = TRUE, 
                            ...) {
  # time.string = gsub(":", "-", gsub(" ", "_", as.character(Sys.time())))
  # if (!is.null(filename)) {
  #   flog.appender(appender.file(paste0(filename, ".", time.string, 
  #                                      ".log")), name = "VennDiagramLogger")
  # }
  # else {
  #   flog.appender(appender.file(paste0("VennDiagram", time.string, 
  #                                      ".log")), name = "VennDiagramLogger")
  # }
  # out.list = as.list(sys.call())
  # out.list[[1]] <- NULL
  # out.string = capture.output(out.list)
  # flog.info(out.string, name = "VennDiagramLogger")
  if (direct.area) {
    if (1 == length(area.vector)) {
      list.names <- category.names
      if (is.null(list.names)) {
        list.names <- ""
      }
      grob.list <- VennDiagram::draw.single.venn(area = area.vector[1], 
                                                 category = list.names, ind = FALSE, ...)
    }
    if (3 == length(area.vector)) {
      grob.list <- VennDiagram::draw.pairwise.venn(area1 = area.vector[1], 
                                                   area2 = area.vector[2], cross.area = area.vector[3], 
                                                   category = category.names, ind = FALSE, print.mode = print.mode, 
                                                   sigdigs = sigdigs, ...)
    }
    if (7 == length(area.vector)) {
      grob.list <- VennDiagram::draw.triple.venn(area1 = 0, 
                                                 area2 = 0, area3 = 0, n12 = 0, n23 = 0, n13 = 0, 
                                                 n123 = 0, category = category.names, ind = FALSE, 
                                                 list.order = 1:3, print.mode = print.mode, sigdigs = sigdigs, 
                                                 area.vector = area.vector, direct.area = TRUE, 
                                                 ...)
    }
    if (15 == length(area.vector)) {
      grob.list <- VennDiagram::draw.quad.venn(area1 = 0, 
                                               area2 = 0, area3 = 0, area4 = 0, n12 = 0, n13 = 0, 
                                               n14 = 0, n23 = 0, n24 = 0, n34 = 0, n123 = 0, 
                                               n124 = 0, n134 = 0, n234 = 0, n1234 = 0, category = category.names, 
                                               ind = FALSE, print.mode = print.mode, sigdigs = sigdigs, 
                                               area.vector = area.vector, direct.area = TRUE, 
                                               ...)
    }
    if (31 == length(area.vector)) {
      grob.list <- VennDiagram::draw.quintuple.venn(area1 = 0, 
                                                    area2 = 0, area3 = 0, area4 = 0, area5 = 0, n12 = 0, 
                                                    n13 = 0, n14 = 0, n15 = 0, n23 = 0, n24 = 0, 
                                                    n25 = 0, n34 = 0, n35 = 0, n45 = 0, n123 = 0, 
                                                    n124 = 0, n125 = 0, n134 = 0, n135 = 0, n145 = 0, 
                                                    n234 = 0, n235 = 0, n245 = 0, n345 = 0, n1234 = 0, 
                                                    n1235 = 0, n1245 = 0, n1345 = 0, n2345 = 0, n12345 = 0, 
                                                    category = category.names, ind = FALSE, print.mode = print.mode, 
                                                    sigdigs = sigdigs, area.vector = area.vector, 
                                                    direct.area = TRUE, ...)
    }
  }
  else {
    if (force.unique) {
      for (i in 1:length(x)) {
        x[[i]] <- unique(x[[i]])
      }
    }
    if ("none" == na) {
      x <- x
    }
    else if ("stop" == na) {
      for (i in 1:length(x)) {
        if (any(is.na(x[[i]]))) {
          flog.error("NAs in dataset", call. = FALSE, 
                     name = "VennDiagramLogger")
          stop("NAs in dataset", call. = FALSE)
        }
      }
    }
    else if ("remove" == na) {
      for (i in 1:length(x)) {
        x[[i]] <- x[[i]][!is.na(x[[i]])]
      }
    }
    else {
      flog.error("Invalid na option: valid options are \"none\", \"stop\", and \"remove\"", 
                 name = "VennDiagramLogger")
      stop("Invalid na option: valid options are \"none\", \"stop\", and \"remove\"")
    }
    if (0 == length(x) | length(x) > 5) {
      flog.error("Incorrect number of elements.", call. = FALSE, 
                 name = "VennDiagramLogger")
      stop("Incorrect number of elements.", call. = FALSE)
    }
    if (1 == length(x)) {
      list.names <- category.names
      if (is.null(list.names)) {
        list.names <- ""
      }
      grob.list <- VennDiagram::draw.single.venn(area = length(x[[1]]), 
                                                 category = list.names, ind = FALSE, ...)
    }
    else if (2 == length(x)) {
      grob.list <- VennDiagram::draw.pairwise.venn(area1 = length(x[[1]]), 
                                                   area2 = length(x[[2]]), cross.area = length(intersect(x[[1]], 
                                                                                                         x[[2]])), category = category.names, ind = FALSE, 
                                                   print.mode = print.mode, sigdigs = sigdigs, ...)
    }
    else if (3 == length(x)) {
      A <- x[[1]]
      B <- x[[2]]
      C <- x[[3]]
      list.names <- category.names
      nab <- intersect(A, B)
      nbc <- intersect(B, C)
      nac <- intersect(A, C)
      nabc <- intersect(nab, C)
      grob.list <- VennDiagram::draw.triple.venn(area1 = length(A), 
                                                 area2 = length(B), area3 = length(C), n12 = length(nab), 
                                                 n23 = length(nbc), n13 = length(nac), n123 = length(nabc), 
                                                 category = list.names, ind = FALSE, list.order = 1:3, 
                                                 print.mode = print.mode, sigdigs = sigdigs, ...)
    }
    else if (4 == length(x)) {
      A <- x[[1]]
      B <- x[[2]]
      C <- x[[3]]
      D <- x[[4]]
      list.names <- category.names
      n12 <- intersect(A, B)
      n13 <- intersect(A, C)
      n14 <- intersect(A, D)
      n23 <- intersect(B, C)
      n24 <- intersect(B, D)
      n34 <- intersect(C, D)
      n123 <- intersect(n12, C)
      n124 <- intersect(n12, D)
      n134 <- intersect(n13, D)
      n234 <- intersect(n23, D)
      n1234 <- intersect(n123, D)
      grob.list <- VennDiagram::draw.quad.venn(area1 = length(A), 
                                               area2 = length(B), area3 = length(C), area4 = length(D), 
                                               n12 = length(n12), n13 = length(n13), n14 = length(n14), 
                                               n23 = length(n23), n24 = length(n24), n34 = length(n34), 
                                               n123 = length(n123), n124 = length(n124), n134 = length(n134), 
                                               n234 = length(n234), n1234 = length(n1234), category = list.names, 
                                               ind = FALSE, print.mode = print.mode, sigdigs = sigdigs, 
                                               ...)
    }
    else if (5 == length(x)) {
      A <- x[[1]]
      B <- x[[2]]
      C <- x[[3]]
      D <- x[[4]]
      E <- x[[5]]
      list.names <- category.names
      n12 <- intersect(A, B)
      n13 <- intersect(A, C)
      n14 <- intersect(A, D)
      n15 <- intersect(A, E)
      n23 <- intersect(B, C)
      n24 <- intersect(B, D)
      n25 <- intersect(B, E)
      n34 <- intersect(C, D)
      n35 <- intersect(C, E)
      n45 <- intersect(D, E)
      n123 <- intersect(n12, C)
      n124 <- intersect(n12, D)
      n125 <- intersect(n12, E)
      n134 <- intersect(n13, D)
      n135 <- intersect(n13, E)
      n145 <- intersect(n14, E)
      n234 <- intersect(n23, D)
      n235 <- intersect(n23, E)
      n245 <- intersect(n24, E)
      n345 <- intersect(n34, E)
      n1234 <- intersect(n123, D)
      n1235 <- intersect(n123, E)
      n1245 <- intersect(n124, E)
      n1345 <- intersect(n134, E)
      n2345 <- intersect(n234, E)
      n12345 <- intersect(n1234, E)
      grob.list <- VennDiagram::draw.quintuple.venn(area1 = length(A), 
                                                    area2 = length(B), area3 = length(C), area4 = length(D), 
                                                    area5 = length(E), n12 = length(n12), n13 = length(n13), 
                                                    n14 = length(n14), n15 = length(n15), n23 = length(n23), 
                                                    n24 = length(n24), n25 = length(n25), n34 = length(n34), 
                                                    n35 = length(n35), n45 = length(n45), n123 = length(n123), 
                                                    n124 = length(n124), n125 = length(n125), n134 = length(n134), 
                                                    n135 = length(n135), n145 = length(n145), n234 = length(n234), 
                                                    n235 = length(n235), n245 = length(n245), n345 = length(n345), 
                                                    n1234 = length(n1234), n1235 = length(n1235), 
                                                    n1245 = length(n1245), n1345 = length(n1345), 
                                                    n2345 = length(n2345), n12345 = length(n12345), 
                                                    category = list.names, ind = FALSE, print.mode = print.mode, 
                                                    sigdigs = sigdigs, ...)
    }
    else {
      flog.error("Invalid size of input object", name = "VennDiagramLogger")
      stop("Invalid size of input object")
    }
  }
  if (length(x) == 2 & !is.null(total.population) & hyper.test) {
    val.p = calculate.overlap.and.pvalue(x[[1]], x[[2]], 
                                         total.population, lower.tail = lower.tail)
    if (is.null(sub)) {
      sub = paste0("p = ", signif(val.p[3], digits = 2))
    }
    else {
      sub = paste0(sub, ", p = ", signif(val.p[3], digits = 2))
    }
  }
  if (!is.null(sub)) {
    grob.list <- add.title(gList = grob.list, x = sub, pos = sub.pos, 
                           fontface = sub.fontface, fontfamily = sub.fontfamily, 
                           col = sub.col, cex = sub.cex)
  }
  if (!is.null(main)) {
    grob.list <- add.title(gList = grob.list, x = main, pos = main.pos, 
                           fontface = main.fontface, fontfamily = main.fontfamily, 
                           col = main.col, cex = main.cex)
  }
  if (!is.null(filename)) {
    current.type <- getOption("bitmapType")
    if (length(grep("Darwin", Sys.info()["sysname"]))) {
      options(bitmapType = "quartz")
    }
    else {
      options(bitmapType = "cairo")
    }
    if ("tiff" == imagetype) {
      tiff(filename = filename, height = height, width = width, 
           units = units, res = resolution, compression = compression)
    }
    else if ("png" == imagetype) {
      png(filename = filename, height = height, width = width, 
          units = units, res = resolution)
    }
    else if ("svg" == imagetype) {
      svg(filename = filename, height = height, width = width)
    }
    else {
      flog.error("You have misspelled your 'imagetype', please try again", 
                 name = "VennDiagramLogger")
      stop("You have misspelled your 'imagetype', please try again")
    }
    grid.draw(grob.list)
    dev.off()
    options(bitmapType = current.type)
    return(1)
  }
  return(grob.list)
}


# run_GO analysis ---------------------------------------------------------

go.asteVector <- function(input){
  output <- c()
  for(i in 1:length(input)){
    output <- paste0(output, ";", input[i])
  }
  output <- substr(x = output, 2, nchar(output))
  output
}

go.intersectGene <- function(A, B){
  intersect(A, B)
}

runTopGO <- function(geneID,GTOGO = NULL, statistic = "fisher", algorithm = "elim",
                     topNodes = 20, database=c("plant","animal"),
                       species = "athaliana_eg_gene", plot = TRUE){
  if(!require(ggplot2)){
    install.packages("ggplot2")
  }
  
  if(!require(biomaRt)){
    BiocManager::install("biomaRt", version = "3.8")
  }
  # 
  if(!require(topGO)){
    BiocManager::install("topGO", version = "3.8")
  }
  if(is.null(GTOGO)){
    if(database == "plant"){
      mart <- useMart(biomart = "plants_mart", dataset = species, host = 'plants.ensembl.org')
    }
    if(database == "animal"){
      mart <- useMart(biomart = "ensembl", dataset = species)
    }
    
    # list all plant GO database
    # listDatasets(mart)
    GTOGO <- getBM(attributes = c( "ensembl_gene_id", "go_id"), mart = mart)
  }
  # GTOGO <- read.table(GO_file,stringsAsFactors = F)
  #GTOGO <- GTOGO[,c(1,6)]
  colnames(GTOGO) <- c( "ensembl_gene_id", "go_id")
  
  GTOGO <- GTOGO[GTOGO$go_id != '', ]
  geneID2GO <- by(GTOGO$go_id, GTOGO$ensembl_gene_id, function(x) as.character(x))
  
  all.genes <- sort(unique(as.character(GTOGO$ensembl_gene_id)))
  int.genes <- geneID
  int.genes <- intersect(int.genes, names(geneID2GO))
  # 1 represent input gene which will be calculate gene enrich
  int.genes <- factor(as.integer(all.genes %in% int.genes))
  names(int.genes) = all.genes
  # 
  go.obj.BP <- new("topGOdata", ontology='BP',
                   allGenes = int.genes,
                   annot = annFUN.gene2GO,
                   gene2GO = geneID2GO)
  
  go.obj.MF <- new("topGOdata", ontology='MF',
                   allGenes = int.genes,
                   annot = annFUN.gene2GO,
                   gene2GO = geneID2GO)
  
  go.obj.CC <- new("topGOdata", ontology='CC',
                   allGenes = int.genes,
                   annot = annFUN.gene2GO,
                   gene2GO = geneID2GO)
  
  ##########retrieve the gene list related to a GO ID######################
  allGO.BP <- genesInTerm(object = go.obj.BP)
  allGO.MF <- genesInTerm(object = go.obj.MF)
  allGO.CC <- genesInTerm(object = go.obj.CC)
  
  #########retrive the significant GO terms
  results.BP <- runTest(go.obj.BP, algorithm = algorithm, statistic = statistic)
  # showSigOfNodes(go.obj.BP, score(results.BP), firstSigNodes = 5, useInfo = 'all')
  results.tab.BP <- GenTable(object = go.obj.BP, elimFisher = results.BP,
                             topNodes = topNodes)
  gene.BP <- genesInTerm(object = go.obj.BP, whichGO = results.tab.BP$GO.ID)
  inter.gene.BP <- lapply(X = gene.BP, FUN = go.intersectGene, B = geneID)
  inter.gene.BP <- unlist(lapply(X = inter.gene.BP, FUN = go.asteVector))
  results.tab.BP$significantGene <- inter.gene.BP
  
  if(length(which(results.tab.BP$elimFisher == "< 1e-30")) != 0){
    results.tab.BP[which(results.tab.BP$elimFisher == "< 1e-30"), ]$elimFisher <- 1e-30
  }
  results.MF <- runTest(go.obj.MF, algorithm = algorithm, statistic = statistic)
  results.tab.MF <- GenTable(object = go.obj.MF, elimFisher = results.MF, 
                             topNodes = topNodes)
  gene.MF <- genesInTerm(object = go.obj.MF, whichGO = results.tab.MF$GO.ID)
  inter.gene.MF <- lapply(X = gene.MF, FUN = go.intersectGene, B = geneID)
  inter.gene.MF <- unlist(lapply(X = inter.gene.MF, FUN = go.asteVector))
  results.tab.MF$significantGene <- inter.gene.MF
  if(length(which(results.tab.MF$elimFisher == "< 1e-30")) != 0){
    results.tab.MF[which(results.tab.MF$elimFisher == "< 1e-30"), ]$elimFisher <- 1e-30
  }
  
  results.CC <- runTest(go.obj.CC, algorithm = algorithm, statistic = statistic)
  results.tab.CC <- GenTable(object = go.obj.CC, elimFisher = results.CC, 
                             topNodes = topNodes)
  gene.CC <- genesInTerm(object = go.obj.CC, whichGO = results.tab.CC$GO.ID)
  inter.gene.CC <- lapply(X = gene.CC, FUN = go.intersectGene, B = geneID)
  inter.gene.CC <- unlist(lapply(X = inter.gene.CC, FUN = go.asteVector))
  results.tab.CC$significantGene <- inter.gene.CC
  if(length(which(results.tab.CC$elimFisher == "< 1e-30")) != 0){
    results.tab.CC[which(results.tab.CC$elimFisher == "< 1e-30"), ]$elimFisher <- 1e-30
  }
  
  
  if(plot){
    df <- data.frame(Category = c(rep("BP", topNodes), rep("CC", topNodes), rep("MF", topNodes)), 
                     x = c(results.tab.BP$Significant, results.tab.CC$Significant, 
                           results.tab.MF$Significant),
                     y = c(-log10(as.numeric(results.tab.BP$elimFisher)), 
                           -log10(as.numeric(results.tab.CC$elimFisher)), 
                           -log10(as.numeric(results.tab.MF$elimFisher))),
                     size = c(-log10(as.numeric(results.tab.BP$elimFisher)),
                              -log10(as.numeric(results.tab.CC$elimFisher)), 
                              -log10(as.numeric(results.tab.MF$elimFisher)))
    )
    
    kk <- ggplot(data = df, aes(x = x, y = y)) + 
      geom_point(aes(color = Category, size = size)) + 
      scale_size_continuous(range = c(2,10)) + 
      labs(x = "The number of significant genes", y = "The adjusted p-values for each GO term")
    print(kk)
  }
  
  results <- list(BP = results.tab.BP, CC = results.tab.CC, MF = results.tab.MF)
  results
}






# KEGG --------------------------------------------------------------------

Ensembl2Entrez <- function(species_name = "Zea mays", geneID = NULL,IDTable=NULL,output_IDtab=F, drop = TRUE){
  if(is.null(IDTable)){
    load("/home/malab10/research/rnaseq_TCtools/analysis_tools/ensembl_species.RData")
    species_type <- species$species[grep(species_name, species$description, ignore.case = TRUE)]
    species_id <- species$dataset[grep(species_name, species$description, ignore.case = TRUE)]
    if(species_type == "ensembl"){
      mart <- useMart(biomart = 'ensembl', dataset = species_id)
    }else{
      biomart <- paste0(species_type, "_mart")
      host <- paste0(species_type, ".ensembl.org")
      mart <- useMart(biomart = biomart, host = host, dataset = species_id)
    }
    IDTable <- getBM(attributes = c("ensembl_gene_id", "entrezgene_id"), mart = mart)
    if(output_IDtab){
      return(IDTable)
    }
  }
  
  idx <- match(geneID, IDTable$ensembl_gene_id)
  entrezID <- IDTable$entrezgene_id[idx]
  
  if(drop){
    entrezID <- unique(na.omit(entrezID))
  }
  entrezID
}
