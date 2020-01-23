analyzePFAM_new <- function(
  switchAnalyzeRlist,
  pathToPFAMresultFile,
  showProgress = TRUE,
  quiet = FALSE,
  pfamMat=NULL) {
  ### Test input
  if(TRUE) {
    if (class(switchAnalyzeRlist) != 'switchAnalyzeRlist') {
      stop(
        'The object supplied to \'switchAnalyzeRlist\' must be a \'switchAnalyzeRlist\''
      )
    }
    if (is.null(switchAnalyzeRlist$orfAnalysis)) {
      stop('ORF needs to be analyzed. Please run analyzeORF and try again.')
    }
    
    
    # qiuzx add
    if (is.null(pfamMat) | ncol(pfamMat) != 16) {
      # file
      if (class(pathToPFAMresultFile) != 'character') {
        stop(
          'The \'pathToPFAMresultFile\' argument must be a string pointing to the PFAM result file(s)'
        )
      }
      if ( ! all(sapply(pathToPFAMresultFile, file.exists)) ) {
        stop('The file(s) \'pathToPFAMresultFile\' points to does not exist')
      }
      
    }else{
      pathToPFAMresultFile <- "none"
    }
    
    
  }
  
  if (showProgress & !quiet) {
    progressBar <- 'text'
  } else {
    progressBar <- 'none'
  }
  
  ### Import result data
  if(TRUE) {
    # qiuzx add 20190919
    if(is.null(pfamMat)){
      ### Test wither headers are included
      temp <-   read.table(
        file = pathToPFAMresultFile[1],
        stringsAsFactors = FALSE,
        fill = TRUE,
        header = FALSE,
        nrows = 1
      )
      if (nrow(temp) == 0) {
        stop('The file pointed to by \'pathToPFAMresultFile\' is empty')
      }
      if (grepl('^<seq|^seq', temp[1, 1])) {
        skipLine <- 1
      } else {
        skipLine <- 0
      }
      
      myPfamResult <- do.call(rbind, plyr::llply(
        pathToPFAMresultFile,
        .fun = function(  aFile  ) {
          read.table(
            file = aFile,
            stringsAsFactors = FALSE,
            fill = TRUE,
            header = FALSE,
            col.names = 1:16,
            skip = skipLine
          )
        }     ))
      
      ### read in pfam resut result
      if (nrow(myPfamResult) == 0) {
        stop('The file pointed to by \'pathToPFAMresultFile\' is empty')
      }
      
      ### Identify which type of PFAM file is supplied
      # Colnames for pfam result
      oldColnames <-
        c(
          'seq_id',
          'alignment_start',
          'alignment_end',
          'envelope_start',
          'envelope_end',
          'hmm_acc',
          'hmm_name',
          'type',
          'hmm_start',
          'hmm_end',
          'hmm_length',
          'bit_score',
          'E_value',
          'significant',
          'clan',
          'residue'
        )
      newColNames <-
        c(
          'seq_id',
          'alignment_start',
          'alignment_end',
          'envelope_start',
          'envelope_end',
          'hmm_acc',
          'hmm_name',
          'hmm_start',
          'hmm_end',
          'hmm_length',
          'bit_score',
          'Individual_E_value',
          'Conditional_E_value',
          'significant',
          'outcompeted',
          'clan'
        )
      
      ### Old style
      if (class(myPfamResult$X8) == 'character') {
        colnames(myPfamResult) <- oldColnames
        myPfamResult$residue[which(myPfamResult$residue == '')] <-NA
        
        ### New style
      } else if (class(myPfamResult$X8) == 'integer') {
        colnames(myPfamResult) <- newColNames
        myPfamResult$clan[which(myPfamResult$clan == '')] <- NA
        
      } else {
        stop('The file(s) supplied is not recogniced as a pfam output.')
      }
      
    }else{
      oldColnames <-
        c(
          'seq_id',
          'alignment_start',
          'alignment_end',
          'envelope_start',
          'envelope_end',
          'hmm_acc',
          'hmm_name',
          'type',
          'hmm_start',
          'hmm_end',
          'hmm_length',
          'bit_score',
          'E_value',
          'significant',
          'clan',
          'residue'
        )
      
      myPfamResult <- pfamMat
      colnames(myPfamResult) <- oldColnames
    }
    
    
    
  }
  
  ### Sanity check that it is a PFAM result file
  if (TRUE) {
    test1 <-      ncol(myPfamResult) == 15 |ncol(myPfamResult) == 16 # the output have 15 or 16 collumns depending on whther active sites are predicted
    test2 <-  all(grepl(      pattern = '^PF|^PB' ,    myPfamResult$hmm_acc,
                              ignore.case = FALSE
    ))                # All pfam hmm starts with PF
    
    if (!all(test1, test2)) {
      stop('The file pointed to by \'pathToPFAMresultFile\' is not a PFAM result file')
    }
    
    # test names
    if (!any(myPfamResult$seq_id %in%
             switchAnalyzeRlist$isoformFeatures$isoform_id)) {
      stop(
        'The transcript ids in the file pointed to by the \'pathToPFAMresultFile\' argument does not match the transcripts stored in the supplied switchAnalyzeRlist'
      )
    }
    
    ### Only add to those with annotated ORF
    #myPfamResult <- myPfamResult[which( myPfamResult$seq_id %in% switchAnalyzeRlist$isoformFeatures$isoform_id ),]
    myPfamResult <- myPfamResult[which(
      myPfamResult$seq_id %in%
        switchAnalyzeRlist$orfAnalysis$isoform_id[which(
          !is.na(switchAnalyzeRlist$orfAnalysis$orfTransciptStart)
        )]
    ), ]
  }
  
  ### Fill in blanks if active residues are included
  if (TRUE) {
    if ('residue' %in% colnames(myPfamResult)) {
      # test whether active residue analysis was performed
      if (any(nchar(myPfamResult$residue) != 0) &
          all(!is.na(myPfamResult$residue))) {
        withActiveRes <- TRUE
      } else {
        withActiveRes <- FALSE
      }
      
      # If active residue analysis was performed do the additional massage of the data
      if (withActiveRes) {
        ### Exchange empty strings in residus with NA's
        myPfamResult$residue[which(!grepl(
          'predicted_active_site',
          myPfamResult$residue
        ))] <- NA
        colnames(myPfamResult)[which(
          colnames(myPfamResult) == 'residue')] <-
          'predicted_active_site'
        
        ### convert active residue to somthing readable
        activeSiteIndex <-
          which(!is.na(myPfamResult$predicted_active_site))
        myPfamResult$predicted_active_site[activeSiteIndex] <-
          sapply(
            myPfamResult$predicted_active_site[activeSiteIndex],
            function(aVec) {
              substr(aVec,
                     start = 23,
                     stop = nchar(aVec) - 1)
            }
          )
      }
    } else {
      withActiveRes <- FALSE
    }
    
    
  }
  
  ### Convert from AA coordinats to transcript and genomic coordinats
  if (TRUE) {
    if (!quiet) {
      message('Converting AA coordinats to transcript and genomic coordinats...')
    }
    ### Remove unwanted columns
    myPfamResult$envelope_start <- NULL
    myPfamResult$envelope_end <- NULL
    myPfamResult$hmm_start <- NULL
    myPfamResult$hmm_end <- NULL
    myPfamResult$hmm_length <- NULL
    
    colnames(myPfamResult)[which(
      grepl('alignment_', colnames(myPfamResult))
    )] <- c('orf_aa_start', 'orf_aa_end')
    
    ### convert from codons to transcript position
    orfStartDF <-
      unique(as.data.frame(
        switchAnalyzeRlist$orfAnalysis[,
                                       c('isoform_id', 'orfTransciptStart')
                                       ]
      ))
    myPfamResult$transcriptStart <-
      (myPfamResult$orf_aa_start  * 3 - 2) +
      orfStartDF[
        match(
          x = myPfamResult$seq_id,
          table = orfStartDF$isoform_id
        ),
        2] - 1
    myPfamResult$transcriptEnd <-
      (myPfamResult$orf_aa_end * 3) +
      orfStartDF[
        match(
          x = myPfamResult$seq_id,
          table = orfStartDF$isoform_id
        ),
        2] - 1
    
    ### convert from transcript to genomic coordinats
    # extract exon data
    myExons <-
      as.data.frame(switchAnalyzeRlist$exons[which(
        switchAnalyzeRlist$exons$isoform_id %in% myPfamResult$seq_id
      ), ])
    myExonsSplit <- split(myExons, f = myExons$isoform_id)
    
    # loop over the individual transcripts and extract the genomic coordiants of the domain and also for the active residues (takes 2 min for 17000 rows)
    myPfamResultDf <-
      plyr::ddply(
        myPfamResult,
        .progress = progressBar,
        .variables = 'seq_id',
        .fun = function(aDF) {
          transcriptId <- aDF$seq_id[1]
          localExons <-
            as.data.frame(myExonsSplit[[transcriptId]])
          
          # extract domain allignement
          localORFalignment <- aDF
          colnames(localORFalignment)[match(
            x = c('transcriptStart', 'transcriptEnd'),
            table = colnames(localORFalignment)
          )] <- c('start', 'end')
          
          # loop over domain alignment (migh be several)
          orfPosList <- list()
          for (j in 1:nrow(localORFalignment)) {
            domainInfo <-
              IsoformSwitchAnalyzeR:::convertCoordinatsTranscriptToGenomic( transcriptCoordinats =  localORFalignment[j, ],
                                                                            exonStructure = localExons              )
            
            ### look into active residues
            if (withActiveRes) {
              if (!is.na(
                localORFalignment$predicted_active_site[j]
              )) {
                activeResInfo <-
                  data.frame(activeRes = as.integer(unlist(
                    strsplit(
                      x = localORFalignment$predicted_active_site[j],
                      split = ','
                    )
                  )))
                
                activeResInfo$start <-
                  activeResInfo$activeRes  * 3 - 2
                activeResInfo$end   <-
                  activeResInfo$activeRes  * 3
                
                activeResInfoList <- list()
                for (k in 1:nrow(activeResInfo)) {
                  activeResInfoList[[as.character(k)]] <-
                    IsoformSwitchAnalyzeR:::convertCoordinatsTranscriptToGenomic(
                      transcriptCoordinats = activeResInfo[k, ],
                      exonStructure = localExons
                    )[, c('pfamStartGenomic',
                          'pfamEndGenomic')]
                }
                activeResInfoDf <-
                  cbind(activeResInfo,
                        do.call(rbind, activeResInfoList))
                
                ### add it to the domain info
                domainInfo$activeResTranscriptStart <-
                  paste(activeResInfoDf$start, collapse = ',')
                domainInfo$activeResTranscriptEnd   <-
                  paste(activeResInfoDf$end, collapse = ',')
                domainInfo$activeResGenomicStart    <-
                  paste(activeResInfoDf$pfamStartGenomic,
                        collapse = ',')
                domainInfo$activeResGenomicEnd      <-
                  paste(activeResInfoDf$pfamEndGenomic,
                        collapse = ',')
              } else {
                ### Add NA instead of residues
                domainInfo$activeResTranscriptStart <- NA
                domainInfo$activeResTranscriptEnd   <- NA
                domainInfo$activeResGenomicStart    <- NA
                domainInfo$activeResGenomicEnd      <- NA
              }
            }
            
            
            orfPosList[[as.character(j)]] <- domainInfo
            
          }
          orfPosDf <- do.call(rbind, orfPosList)
          
          return(cbind(aDF, orfPosDf))
        }
      )
    
  }
  
  ### Add analysis to switchAnalyzeRlist
  if (TRUE) {
    ### reorder data.frame
    # make sure the basic data is last
    colnames(myPfamResultDf)[1] <- 'isoform_id'
    newOrderNames <- c('isoform_id', 'hmm_acc', 'hmm_name', 'clan')
    myPfamResultDf <-
      myPfamResultDf[, c(
        which(colnames(myPfamResultDf) %in% newOrderNames) ,
        which(!colnames(myPfamResultDf) %in% newOrderNames)
      )]
    
    # sort
    myPfamResultDf <-
      myPfamResultDf[order(
        myPfamResultDf$isoform_id,
        myPfamResultDf$transcriptStart,
        myPfamResultDf$hmm_name
      ), ]
    
    # if active residues are pressent put them last
    if (withActiveRes) {
      newOrder <-
        c(which(
          !grepl(
            'Predicted_active_site|ActiveRes',
            colnames(myPfamResultDf)
          )
        ), which(
          grepl(
            'Predicted_active_site|ActiveRes',
            colnames(myPfamResultDf)
          )
        ))
      myPfamResultDf <- myPfamResultDf[, newOrder]
    }
    
    #myPfamResultDf$pfamStarExon <- NULL
    #myPfamResultDf$pfamEndExon <- NULL
    
    # add the pfam results to the switchAnalyzeRlist object
    switchAnalyzeRlist$domainAnalysis <- myPfamResultDf
    
    # add indication to transcriptDf
    switchAnalyzeRlist$isoformFeatures$domain_identified <- 'no'
    switchAnalyzeRlist$isoformFeatures$domain_identified[which(
      is.na(switchAnalyzeRlist$isoformFeatures$PTC)
    )] <- NA # sets NA for those not analyzed
    
    switchAnalyzeRlist$isoformFeatures$domain_identified[which(
      switchAnalyzeRlist$isoformFeatures$isoform_id %in%
        myPfamResultDf$isoform_id
    )] <- 'yes'
  }
  
  n <- length(unique(myPfamResultDf$isoform_id))
  p <-
    round(n / length(unique(
      switchAnalyzeRlist$isoformFeatures$isoform_id
    )) * 100, digits = 2)
  
  if (!quiet) {
    message(paste(
      'Added domain information to ',
      n,
      ' (',
      p,
      '%) transcripts',
      sep = ''
    ))
  }
  return(switchAnalyzeRlist)
}

importCufflinksFiles_new <- function(
  pathToGTF,
  pathToGeneDEanalysis,
  pathToIsoformDEanalysis,
  pathToGeneFPKMtracking,
  pathToIsoformFPKMtracking,
  pathToIsoformReadGroupTracking,
  pathToSplicingAnalysis = NULL,
  pathToReadGroups,
  pathToRunInfo,
  fixCufflinksAnnotationProblem = TRUE,
  addIFmatrix = TRUE,
  quiet = FALSE
) {
  ### Test that files exist
  if (TRUE) {
    if( pathToGTF == '' ) {
      stop(
        paste(
          'The \'pathToGTF\' argument does not lead anywhere (acutally you just suppled "" to the argument).',
          '\nDid you try to use the system.file("your/quant/dir/", package="IsoformSwitchAnalyzeR")',
          'to import your own data? The system.file() should only be used',
          'to access the example data stored in the IsoformSwitchAnalyzeR package.',
          'To access your own data simply provide the string to the directory with the data as:',
          '"path/to/quantification/".',
          sep=' '
        )
      )
    }
    if( pathToGeneDEanalysis == '' ) {
      stop(
        paste(
          'The \'pathToGeneDEanalysis\' argument does not lead anywhere (acutally you just suppled "" to the argument).',
          '\nDid you try to use the system.file("your/quant/dir/", package="IsoformSwitchAnalyzeR")',
          'to import your own data? The system.file() should only be used',
          'to access the example data stored in the IsoformSwitchAnalyzeR package.',
          'To access your own data simply provide the string to the directory with the data as:',
          '"path/to/quantification/".',
          sep=' '
        )
      )
    }
    
    # pathToGTF
    if (!file.exists(pathToGTF)) {
      stop('The \'pathToGTF\' argument does not point to an acutal file')
    }
    # DE
    if (!file.exists(pathToGeneDEanalysis))    {
      stop('The \'pathToGeneDEanalysis\' argument does not point to an acutal file')
    }
    if (!file.exists(pathToIsoformDEanalysis)) {
      stop('The \'pathToIsoformDEanalysis\' argument does not point to an acutal file')
    }
    # Tracking
    if (!file.exists(pathToGeneFPKMtracking))    {
      stop('The \'pathToGeneFPKMtracking\' argument does not point to an acutal file')
    }
    if (!file.exists(pathToIsoformFPKMtracking)) {
      stop(
        'The \'pathToIsoformFPKMtracking\' argument does not point to an acutal file'
      )
    }
    if (!file.exists(pathToIsoformReadGroupTracking)) {
      stop(
        'The \'pathToIsoformReadGroupTracking\' argument does not point to an acutal file'
      )
    }
    # splicing
    if (!is.null(pathToSplicingAnalysis)) {
      if (!file.exists(pathToSplicingAnalysis)) {
        stop(
          'The \'pathToSplicingAnalysis\' argument does not point to an acutal file'
        )
      }
    }
    # info
    if (!file.exists(pathToReadGroups)) {
      stop('The \'pathToReadGroups\' argument does not point to an acutal file')
    }
    if (!file.exists(pathToRunInfo))    {
      stop('The \'pathToRunInfo\' argument does not point to an acutal file')
    }
  }
  
  ### Import the supplied files (not gtf)
  if (TRUE) {
    if (!quiet) { message('Step 1 of 5: Importing data...')}
    suppressMessages(
      geneDiffanalysis     <-
        readr::read_tsv(
          file = pathToGeneDEanalysis,
          col_names = TRUE
        )
    )
    suppressMessages(
      isoformDiffanalysis  <-
        readr::read_tsv(
          file = pathToIsoformDEanalysis,
          col_names = TRUE
        )
    )
    suppressMessages(
      geneAnnotation       <-
        readr::read_tsv(
          file = pathToGeneFPKMtracking,
          col_names = TRUE
        )
    )
    suppressMessages(
      isoformAnnotation    <-
        readr::read_tsv(
          file = pathToIsoformFPKMtracking,
          col_names = TRUE
        )
    )
    suppressMessages(
      isoRepExp        <-
        read.table(
          file = pathToIsoformReadGroupTracking,
          header = TRUE,
          sep='\t',
          stringsAsFactors = FALSE
        )
    )
    suppressMessages(
      cuffSplicing         <-
        readr::read_tsv(
          file = pathToSplicingAnalysis,
          col_names = TRUE
        )
    )
    
    suppressMessages(
      readGroup <-
        read.table(
          file = pathToReadGroups,
          sep='\t',
          header = TRUE
        )
    )
    suppressMessages(
      runInfo   <-
        readr::read_tsv(
          file = pathToRunInfo,
          col_names = TRUE
        )
    )
  }
  
  ### "Test" that the data.files are what they are supposed to be
  if (TRUE) {
    ### gene diff analysis
    q1 <-
      !all(
        colnames(geneDiffanalysis) %in% c(
          "test_id",
          "gene_id",
          "gene",
          "locus",
          "sample_1",
          "sample_2",
          "status",
          "value_1",
          "value_2",
          "log2(fold_change)",
          "test_stat",
          "p_value",
          "q_value",
          "significant"
        )
      )
    if (q1) {
      stop(paste(
        'The file supplied to pathToGeneDEanalysis does not appear',
        'to be the result of the CuffDiff gene expression analysis.'
      ))
    }
    ### transcript diff analysis
    q1 <-
      !all(
        colnames(isoformDiffanalysis) %in% c(
          "test_id",
          "gene_id",
          "gene",
          "locus",
          "sample_1",
          "sample_2",
          "status",
          "value_1",
          "value_2",
          "log2(fold_change)",
          "test_stat",
          "p_value",
          "q_value",
          "significant"
        )
      )
    if (q1) {
      stop(paste(
        'The file supplied to isoformDiffanalysis does not appear to',
        'be the result of the CuffDiff transcript expression analysis.'
      ))
    }
    
    q2 <-
      sum(grepl(
        'TCONS', isoformDiffanalysis$test_id
      )) != nrow(isoformDiffanalysis)
    if (q2) {
      warning(paste(
        'It looks like you have NOT been doing transcript\n',
        'reconstruction/assembly with Cufflinks/Cuffdiff.\n',
        'If you have not reconstructed transcripts we receomend to use Kallisto or Salmon\n',
        'to do the quantification instead - they are more accurate and have better biase correction methods.'
      ))
    }
    
    
    
    ### gene annoation
    q1 <-
      !all(
        colnames(geneAnnotation)[1:8] %in% c(
          "tracking_id",
          "class_code",
          "nearest_ref_id",
          "gene_id",
          "gene_short_name",
          "tss_id",
          "locus",
          "length"
        )
      )
    if (q1) {
      stop(paste(
        'The file supplied to geneAnnotation does not appear to be the',
        'gene FPKM traccking of the CuffDiff gene FPKM trascking analysis.'
      ))
    }
    ### transcript annoation
    q1 <-
      !all(
        colnames(isoformAnnotation)[1:8] %in% c(
          "tracking_id",
          "class_code",
          "nearest_ref_id",
          "gene_id",
          "gene_short_name",
          "tss_id",
          "locus",
          "length"
        )
      )
    if (q1) {
      stop(paste(
        'The file supplied to isoformAnnotation does not appear to be',
        'the isoform FPKM tracking of the CuffDiff transcript analysis.'
      ))
    }
    
    ### rep expression
    q1 <-
      !all(
        colnames(isoRepExp)[1:4] %in% c(
          "tracking_id", "condition", "replicate", "raw_frags"
        )
      )
    if (q1) {
      stop(paste(
        'The file supplied to pathToIsoformCountTracking does not',
        'appear to be the isoform count tracking of the CuffDiff',
        'transcript analysis.'
      ))
    }
    
    ### splicing analysis
    q1 <-
      !all(
        colnames(cuffSplicing) %in% c(
          "test_id",
          "gene_id",
          "gene",
          "locus",
          "sample_1",
          "sample_2",
          "status",
          "value_1",
          "value_2",
          "sqrt(JS)",
          "test_stat",
          "p_value",
          "q_value",
          "significant"
        )
      )
    if (q1) {
      stop(
        'The file supplied to cuffSplicing does not appear to be the',
        'result of the CuffDiff differential analysis of alternative splicing.'
      )
    }
    
    ### Read grous
    q1 <-
      !all(
        colnames(readGroup) %in% c(
          "file",
          "condition",
          "replicate_num",
          "total_mass",
          "norm_mass",
          "internal_scale",
          "external_scale"
        )
      )
    q2 <-
      !all(readGroup$condition %in% unique(
        unlist(geneDiffanalysis[, c('sample_1', 'sample_2')]))
      )
    if (q1 | q2) {
      stop(paste(
        'The file supplied to readGroup does not appear to be the',
        'pathToReadGroups of the CuffDiff transcript analysis.'
      ))
    }
    
    ### Run info
    q1 <- !all(colnames(runInfo) %in% c("param", 'value'))
    q2 <-
      !all(
        runInfo$param %in% c(
          "cmd_line",
          "version",
          "SVN_revision",
          "boost_version"
        )
      )
    if (q1 | q2) {
      stop(paste(
        'The file supplied to runInfo does not appear to be',
        'the runInfo of the CuffDiff transcript analysis.'
      ))
    }
    
  }
  
  ### Massage and merge gene and isoform annoation and DE analysis
  if (TRUE) {
    if (!quiet) { message('Step 2 of 5: Merging gene and isoform expression...')}
    ### Design matrix
    readGroup$sample_name <-  stringr::str_c(readGroup$condition, '_', readGroup$replicate_num)
    designMatrix <- readGroup[, c('sample_name', 'condition')]
    colnames(designMatrix) <- c('sampleID', 'condition')
    
    ### Massage data frames
    if (TRUE) {
      # gene
      geneDiffanalysis <-    geneDiffanalysis[, -which(     colnames(geneDiffanalysis) %in%  c('test_id', 'gene', 'locus', 'test_stat')    )]
      colnames(geneDiffanalysis)[4:ncol(geneDiffanalysis)] <-
        paste(
          "gene_",
          colnames(geneDiffanalysis)[4:ncol(geneDiffanalysis)] ,
          sep = "") # add gene to the colnames so they can be destinquished from the gene diff data
      colnames(geneDiffanalysis) <- gsub(
        'gene_log2.fold_change.',
        'gene_log2_fold_change',
        colnames(geneDiffanalysis)
      )
      
      # info
      colnames(isoformAnnotation)[1] <- 'isoform_id'
      isoformAnnotation2 <-
        isoformAnnotation[, na.omit(match(
          c(
            'isoform_id',
            'gene_id',
            'gene_short_name',
            'nearest_ref_id',
            'class_code',
            'tss_id',
            'CDS_id',
            'length',
            'locus'
          ),
          colnames(isoformAnnotation)
        ))]
      
      colnames(isoformAnnotation2)[which(
        colnames(isoformAnnotation2) == 'gene_short_name'
      )] <- 'gene_name'
      
      # iso
      isoformDiffanalysis <- isoformDiffanalysis[, which(
        !colnames(isoformDiffanalysis) %in%
          c('gene_id', 'gene', 'test_stat')
      )]
      colnames(isoformDiffanalysis)[5:ncol(isoformDiffanalysis)] 	<-
        paste(
          "iso_",
          colnames(isoformDiffanalysis)[5:ncol(isoformDiffanalysis)],
          sep = "") # add gene to the colnames so they can be destinquished from the gene diff data
      colnames(isoformDiffanalysis)[1] <- 'isoform_id'
      colnames(isoformDiffanalysis) <-
        gsub(
          'iso_log2.fold_change.',
          'iso_log2_fold_change',
          colnames(isoformDiffanalysis)
        )
      
      # rep expression
      isoRepExp2 <-
        isoRepExp[, c("tracking_id",
                      "condition",
                      "replicate",
                      "raw_frags")]
      colnames(isoRepExp2)[1] <- 'isoform_id'
      isoRepExp2$rep <-       paste(isoRepExp2$condition, isoRepExp2$replicate, sep = '_')
      isoRepExp2 <-    reshape2::dcast(data = isoRepExp2,   isoform_id ~ rep,  value.var = 'raw_frags')
      
      
      ### rep fpkm
      # iso
      isoRepFpkm <- isoRepExp[, c(
        "tracking_id",
        "condition",
        "replicate",
        "FPKM"
      )]
      colnames(isoRepFpkm)[1] <- 'isoform_id'
      isoRepFpkm$rep <-
        paste(isoRepFpkm$condition, isoRepFpkm$replicate, sep = '_')
      isoRepFpkm <-
        reshape2::dcast(data = isoRepFpkm,
                        isoform_id ~ rep,
                        value.var = 'FPKM')
      
      ### Gene
      isoRepFpkm2 <- isoRepFpkm
      isoRepFpkm2$gene_id <- isoformAnnotation2$gene_id[match(  isoRepFpkm2$isoform_id, isoformAnnotation2$isoform_id      )]
      isoRepFpkm2$isoform_id <- NULL
      
      geneRepFpkm <- isoformToGeneExp(isoRepFpkm2, quiet = TRUE)
      
      ### Calculate means
      rownames(isoRepFpkm) <- isoRepFpkm$isoform_id
      isoRepFpkm$isoform_id <- NULL
      isoMean <- rowMeans(isoRepFpkm)
      
      # qiuzx delete
      # rownames(geneRepFpkm) <- geneRepFpkm$gene_id
      # geneRepFpkm$gene_id <- NULL
      geneMean <- rowMeans(geneRepFpkm)
      
      ### add means
      geneDiffanalysis$gene_overall_mean <- geneMean[match(     geneDiffanalysis$gene_id, names(geneMean)  )]
      
      isoformDiffanalysis$iso_overall_mean <- isoMean[match(
        isoformDiffanalysis$isoform_id, names(isoMean)
      )]
      
    }
    
    ### Extract standard error
    if (TRUE) {
      ### Tjek if the Isoform CI collums are switches
      ciLowColumn <-
        which(grepl('_conf_lo', colnames(isoformAnnotation)))[1]
      ciHighColumn <-
        which(grepl('_conf_hi', colnames(isoformAnnotation)))[1]
      
      if (all(isoformAnnotation[, ciHighColumn] >= isoformAnnotation[, ciLowColumn])) {
        highString <- '_conf_hi'
      } else {
        highString <- '_conf_lo'
      }
      
      ### extract isoform sddev from CI
      # fpkm
      isoformFPKM <-
        isoformAnnotation[, which(grepl(
          'isoform_id|_FPKM', colnames(isoformAnnotation)
        ))]
      isoformFPKM <- reshape2::melt(isoformFPKM, id.vars = 'isoform_id')
      isoformFPKM$variable <-
        gsub('_FPKM$', '', isoformFPKM$variable)
      colnames(isoformFPKM)[3] <- 'expression'
      # ci high
      isoformFPKMciHi <-
        isoformAnnotation[, which(grepl(
          paste('isoform_id|', highString, sep = ''),
          colnames(isoformAnnotation)
        ))]
      isoformFPKMciHi <-
        reshape2::melt(isoformFPKMciHi, id.vars = 'isoform_id')
      isoformFPKMciHi$variable <-
        gsub(highString, '', isoformFPKMciHi$variable)
      colnames(isoformFPKMciHi)[3] <- 'ci_hi'
      # stderr
      isoformFPKMcombined <-
        dplyr::inner_join(isoformFPKM,
                          isoformFPKMciHi,
                          by = c('isoform_id', 'variable'))
      isoformFPKMcombined$iso_stderr <-
        (isoformFPKMcombined$ci_hi - isoformFPKMcombined$expression) / 2 # How it's done in cufflinks source code
      isoformFPKMcombined$expression <- NULL
      isoformFPKMcombined$ci_hi <- NULL
      colnames(isoformFPKMcombined) <-
        c('isoform_id', 'sample_name', 'iso_stderr')
      
      ### Tjek if the gene CI collums are switches
      ciLowColumn <-
        which(grepl('_conf_lo', colnames(geneAnnotation)))[1]
      ciHighColumn <-
        which(grepl('_conf_hi', colnames(geneAnnotation)))[1]
      
      if (all(
        geneAnnotation[, ciHighColumn] >= geneAnnotation[, ciLowColumn])
      ) {
        highString <- '_conf_hi'
      } else {
        highString <- '_conf_lo'
      }
      
      ### extract gene sddev from CI
      # fpkm
      geneFPKM <- geneAnnotation[, which(grepl(
        'tracking_id|_FPKM', colnames(geneAnnotation)
      ))]
      geneFPKM <- reshape2::melt(geneFPKM, id.vars = 'tracking_id')
      geneFPKM$variable <-
        gsub('_FPKM$', '', geneFPKM$variable)
      colnames(geneFPKM)[3] <- 'expression'
      # ci high
      geneFPKMciHi <- geneAnnotation[, which(grepl(
        paste('tracking_id|', highString, sep = ''),
        colnames(geneAnnotation)
      ))]
      geneFPKMciHi <- reshape2::melt(geneFPKMciHi, id.vars = 'tracking_id')
      geneFPKMciHi$variable <- gsub(highString, '', geneFPKMciHi$variable)
      colnames(geneFPKMciHi)[3] <- 'ci_hi'
      # stderr
      geneFPKMcombined <-
        dplyr::inner_join(geneFPKM,
                          geneFPKMciHi,
                          by = c('tracking_id', 'variable'))
      geneFPKMcombined$iso_stderr <-
        (geneFPKMcombined$ci_hi - geneFPKMcombined$expression) / 2 # how it's done in cufflinks sourece code
      geneFPKMcombined$expression <- NULL
      geneFPKMcombined$ci_hi <- NULL
      colnames(geneFPKMcombined) <-
        c('gene_id', 'sample_name', 'gene_stderr')
      
      
      ## Merge stderr with DE analysis
      #isoformDiffanalysis <-
      #    merge(
      #        isoformDiffanalysis,
      #        isoformFPKMcombined,
      #        by.x = c('isoform_id', 'sample_2'),
      #        by.y = c('isoform_id', 'sample_name')
      #    )
      isoformDiffanalysis <- suppressWarnings( dplyr::inner_join(
        isoformDiffanalysis,
        isoformFPKMcombined,
        by=c("sample_2" = "sample_name", "isoform_id" = "isoform_id")
      ) )
      colnames(isoformDiffanalysis)[which( grepl(
        'iso_stderr', colnames(isoformDiffanalysis))
      )] <- 'iso_stderr_2'
      
      #isoformDiffanalysis <- merge(
      #        isoformDiffanalysis,
      #        isoformFPKMcombined,
      #        by.x = c('isoform_id', 'sample_1'),
      #        by.y = c('isoform_id', 'sample_name')
      #    )
      isoformDiffanalysis <- suppressWarnings( dplyr::inner_join(
        isoformDiffanalysis,
        isoformFPKMcombined,
        by=c("sample_2" = "sample_name", "isoform_id" = "isoform_id")
      ) )
      colnames(isoformDiffanalysis)[which(grepl(
        'iso_stderr$',
        colnames(isoformDiffanalysis),
        perl = TRUE
      ))] <- c('iso_stderr_1')
      
      isoformDiffanalysis <-
        isoformDiffanalysis[, c(
          'isoform_id',
          'sample_1',
          'sample_2',
          'iso_status',
          'iso_overall_mean',
          'iso_value_1',
          'iso_value_2',
          'iso_stderr_1',
          'iso_stderr_2',
          'iso_log2_fold_change',
          'iso_p_value',
          'iso_q_value',
          'iso_significant'
        )]
      
      ### Extract and add gene stderr
      #geneDiffanalysis <-
      #    merge(
      #        geneDiffanalysis,
      #        geneFPKMcombined,
      #        by.x = c('gene_id', 'sample_2'),
      #        by.y = c('gene_id', 'sample_name')
      #    )
      geneDiffanalysis <- suppressWarnings( dplyr::inner_join(
        geneDiffanalysis,
        geneFPKMcombined,
        by=c("sample_2" = "sample_name", "gene_id" = "gene_id")
      ) )
      colnames(geneDiffanalysis)[ which(grepl(
        'gene_stderr', colnames(geneDiffanalysis)
      ))] <- 'gene_stderr_2'
      
      #geneDiffanalysis <- merge(
      #        geneDiffanalysis,
      #        geneFPKMcombined,
      #        by.x = c('gene_id', 'sample_1'),
      #        by.y = c('gene_id', 'sample_name')
      #    )
      geneDiffanalysis <- suppressWarnings( dplyr::inner_join(
        geneDiffanalysis,
        geneFPKMcombined,
        by=c("sample_1" = "sample_name", "gene_id" = "gene_id")
      ) )
      colnames(geneDiffanalysis)[which(grepl(
        'gene_stderr$',
        colnames(geneDiffanalysis),
        perl = TRUE
      ))] <- c('gene_stderr_1')
      
      geneDiffanalysis <-
        geneDiffanalysis[, c(
          'gene_id',
          'sample_1',
          'sample_2',
          'gene_status',
          'gene_overall_mean',
          'gene_value_1',
          'gene_value_2',
          'gene_stderr_1',
          'gene_stderr_2',
          'gene_log2_fold_change',
          'gene_p_value',
          'gene_q_value',
          'gene_significant'
        )]
      
    }
    
    
    ### Merge data
    if (TRUE) {
      ### Meger gene DE and annotation data
      isoformData <-
        dplyr::inner_join(isoformAnnotation2, geneDiffanalysis, by = 'gene_id')
      
      ### Merge with iso DE
      isoformData <-
        dplyr::inner_join(
          isoformData,
          isoformDiffanalysis,
          by = c('isoform_id', 'sample_1', 'sample_2')
        )
      
      ### Massage again
      colnames(isoformData)[which(
        colnames(isoformData) == 'tss_id'
      )] <- 'TSS_group_id'
      
    }
    
  }
  
  ### Obtain transcript structure information
  if (TRUE) {
    if (!quiet) { message('Step 3 of 5: Obtaining annotation...')}
    
    ### Import file
    exonFeatures <-  rtracklayer::import(pathToGTF, format = 'gtf')
    if (length(exonFeatures) == 0)
      stop("No exon information extracted from GTF")
    
    ### Filter for what is needed
    exonFeatures <- exonFeatures[
      which(tolower(exonFeatures$type) == 'exon'),
      c('gene_id', 'transcript_id')
      ]
    
    ### rename
    colnames(exonFeatures@elementMetadata) <- gsub(
      'transcript_id', 'isoform_id',
      colnames(exonFeatures@elementMetadata)
    )
  }
  
  ### Check it is the same transcripts in transcript structure and expression info
  if (TRUE) {
    myUnion     <-
      unique(c(
        isoformData$isoform_id,
        exonFeatures$isoform_id,
        isoRepExp2$isoform_id
      ))
    myIntersect <- intersect(
      intersect(isoformData$isoform_id, exonFeatures$isoform_id),
      isoRepExp2$isoform_id
    )
    
    # If there are descripencies
    if(length(myIntersect) == 0) {
      stop(
        paste(
          'No overlap between isoform annotation',
          'and isoform expression data was found',
          sep=' '
        )
      )
    }
    
    
    if (length(myUnion) != length(myIntersect)) {
      isoformData <- isoformData[which(
        isoformData$isoform_id     %in% myIntersect), ]
      exonFeatures <- exonFeatures[which(
        exonFeatures$isoform_id    %in% myIntersect), ]
      isoRepExp2 <- isoRepExp2[which(
        isoRepExp2$isoform_id       %in% myIntersect), ]
      
      if (!quiet) {
        message(
          paste(
            'There were discrepencies between the GTF and the',
            'expression analysis files. To solve this',
            abs(length(myUnion) - length(myIntersect)) ,
            'transcripts were removed.',
            sep = ' '
          )
        )
      }
    }
  }
  
  ### Fix to correct for Cufflinks annotation problem where cufflinks assignes
  # transcripts from several annotated genes to 1 cuffgene
  if (   fixCufflinksAnnotationProblem ) {
    if (!quiet) { message('Step 4 of 3: Fixing cufflinks annotation problem...')}
    
    geneName <- unique(isoformData[, c('gene_id', 'gene_name')])
    geneNameSplit <-
      split(geneName$gene_name , f = geneName$gene_id)
    # remove all unique
    geneNameSplit <-
      geneNameSplit[which(sapply(geneNameSplit, function(x)
        length(unique(x))) > 1)]
    
    if (length(geneNameSplit) > 0) {
      # if there are any problems
      #get indexes of those affected
      geneNameIndexesData     <-
        which(isoformData$gene_id %in% names(geneNameSplit))
      geneNameIndexesFeatures <-
        which(exonFeatures$spliceR.gene_id %in% names(geneNameSplit))
      
      # combine names of cuffgenes and gene short name
      isoformData$gene_id[geneNameIndexesData]              <-
        paste(isoformData$gene_id[geneNameIndexesData]        ,
              isoformData$gene_name[geneNameIndexesData],
              sep = ':')
      exonFeatures$spliceR.gene_id[geneNameIndexesFeatures] <-
        paste(exonFeatures$spliceR.gene_id[geneNameIndexesFeatures],
              exonFeatures$spliceR.gene_name[geneNameIndexesFeatures],
              sep = ':')
      
      ## Correct gene expression levels and differntial analysis
      problematicGenes <-
        isoformData[geneNameIndexesData, c(
          'isoform_id',
          'gene_id',
          'sample_1',
          'sample_2',
          'gene_overall_mean',
          'gene_value_1',
          'gene_value_2',
          'gene_stderr_1',
          'gene_stderr_2',
          'gene_log2_fold_change',
          'gene_p_value',
          'gene_q_value',
          'gene_significant',
          'iso_status',
          'iso_overall_mean',
          'iso_value_1',
          'iso_value_2'
        )]
      problematicGenesSplit <-
        split(problematicGenes, f = problematicGenes[
          ,c('gene_id', 'sample_1', 'sample_2')], drop =TRUE)
      
      correctedGenes <-
        plyr::ldply(
          problematicGenesSplit,
          .fun = function(df) {
            # df <- problematicGenesSplit[[1]]
            df$gene_overall_mean <- sum(df$iso_overall_mean)
            df$gene_value_1 <- sum(df$iso_value_1)
            df$gene_value_2 <- sum(df$iso_value_2)
            df$gene_stderr_1 <- NA
            df$gene_stderr_2 <- NA
            df$gene_log2_fold_change <- log2(
              (df$gene_value_2[2] + 0.0001) /
                (df$gene_value_1[1] + 0.0001)
            )
            df$gene_p_value <- 1
            df$gene_q_value <- 1
            df$gene_significant <- 'no'
            df$iso_status <- 'NOTEST'
            return(df)
          }
        )
      
      # sort so genes end up being in correct order for overwriting
      correctedGenes <-
        correctedGenes[order(
          correctedGenes$isoform_id,
          correctedGenes$gene_id,
          correctedGenes$sample_1,
          correctedGenes$sample_2
        ), -1] # -1 removes the index created by ldply
      # overwrite problematic genes
      isoformData[geneNameIndexesData, c(
        'gene_id',
        'sample_1',
        'sample_2',
        'gene_overall_mean',
        'gene_value_1',
        'gene_value_2',
        'gene_stderr_1',
        'gene_stderr_2',
        'gene_log2_fold_change',
        'gene_p_value',
        'gene_q_value',
        'gene_significant',
        'iso_status',
        'iso_overall_mean',
        'iso_value_1',
        'iso_value_2'
      )] <- correctedGenes[, -1] # -1 removes the isoform id
      
      
      ### Add to exons
      exonFeatures$gene_id <- isoformData$gene_id[match(exonFeatures$isoform_id, isoformData$isoform_id   )]
      
      if (!quiet) {
        message(
          paste(
            "    Cufflinks annotation problem was fixed for",
            length(geneNameSplit),
            "Cuff_genes",
            sep = ' '
          )
        )
      }
    } else {
      if (!quiet) {
        message(paste(
          'No instances of a Cufflinks annotation',
          'problem found - no changes were made'
        ))
      }
    }
  } # end of fix cufflinks annotatopn problem
  if ( ! fixCufflinksAnnotationProblem ) {
    if (!quiet) { message('Step 4 of 5: Skipped fixing of cufflinks annotation problem (due to fixCufflinksAnnotationProblem argument)...')}
  }
  
  if (!quiet) { message('Step 5 of 5: Creating switchAanalyzeRlist...')}
  
  ### Calculate IF values
  localAnnot <- unique(isoformData[,c('gene_id','isoform_id')])
  ## qiuzx add
  com_isoID <- intersect(rownames(isoRepFpkm),localAnnot$isoform_id)
  isoRepFpkm <- isoRepFpkm[com_isoID,]
  localAnnot <- localAnnot[localAnnot$isoform_id %in% com_isoID,]
  ##
  ifMatrix <- isoformToIsoformFraction(
    isoformRepExpression = isoRepFpkm,
    isoformGeneAnnotation = localAnnot,
    quiet = TRUE
  )
  
  ### Summarize IF
  myMeanIF <- rowMeans(ifMatrix[,designMatrix$sampleID,drop=FALSE], na.rm = TRUE)
  ifMeanDf <- plyr::ddply(
    .data = designMatrix,
    .variables = 'condition',
    .fun = function(aDF) { # aDF <- switchAnalyzeRlist$designMatrix[1:2,]
      tmp <- rowMeans(ifMatrix[,aDF$sampleID,drop=FALSE], na.rm = TRUE)
      data.frame(
        isoform_id=names(tmp),
        mean=tmp,
        stringsAsFactors = FALSE
      )
    }
  )
  
  isoformData$IF_overall <- myMeanIF[match(
    isoformData$isoform_id, names(myMeanIF)
  )]
  isoformData$IF1 <- ifMeanDf$mean[match(
    stringr::str_c(isoformData$isoform_id,isoformData$sample_1),
    stringr::str_c(ifMeanDf$isoform_id, ifMeanDf$condition)
  )]
  isoformData$IF2 <- ifMeanDf$mean[match(
    stringr::str_c(isoformData$isoform_id,isoformData$sample_2),
    stringr::str_c(ifMeanDf$isoform_id, ifMeanDf$condition)
  )]
  isoformData$dIF <- isoformData$IF2 - isoformData$IF1
  
  ### Add q-values
  if (!is.null(pathToSplicingAnalysis)) {
    ### Add cufflinks analysis isoform switch analysis results
    isoformData$isoform_switch_q_value <- NA
    isoformData$gene_switch_q_value <-
      cuffSplicing$q_value[match(
        paste(
          isoformData$gene_id,
          isoformData$sample_1,
          isoformData$sample_2,
          sep = '_'
        ),
        paste(
          cuffSplicing$gene_id,
          cuffSplicing$sample_1,
          cuffSplicing$sample_2,
          sep = '_'
        )
      )]
  } else {
    ### Add collumns for isoform switch analysis results
    isoformData$isoform_switch_q_value <- NA
    isoformData$gene_switch_q_value <- NA
  }
  
  ### Reorder a bit
  ofInterest <- c('isoform_id','gene_id','gene_name','sample_1','sample_2')
  isoformData <- isoformData[, c(
    match( ofInterest, colnames(isoformData)),
    which( ! colnames(isoformData) %in% ofInterest)
  )]
  colnames(isoformData)[4:5] <- c('condition_1', 'condition_2')
  isoformData <- as.data.frame(isoformData)
  
  ### Extract run info
  # cufflinks version
  cuffVersion <- runInfo$value[2]
  
  # replicate numbers
  nrRep <- table(readGroup$condition)
  nrRep <-
    data.frame(
      condition = names(nrRep),
      nrReplicates = as.vector(nrRep),
      row.names = NULL,
      stringsAsFactors = FALSE
    )
  
  
  isoRepFpkm$isoform_id <- rownames(isoRepFpkm)
  rownames(isoRepFpkm) <- NULL
  isoRepFpkm <- isoRepFpkm[,c(
    which(colnames(isoRepFpkm) == 'isoform_id'),
    which(colnames(isoRepFpkm) != 'isoform_id')
  )]
  
  # Return SpliceRList
  switchAnalyzeRlist <- createSwitchAnalyzeRlist(
    isoformFeatures = isoformData,
    exons = exonFeatures,
    designMatrix = designMatrix,
    isoformCountMatrix = isoRepExp2,
    isoformRepExpression = isoRepFpkm,
    sourceId = paste("cufflinks", cuffVersion , sep = '_')
  )
  
  if (!is.null(pathToSplicingAnalysis) & nrow(cuffSplicing)) {
    switchAnalyzeRlist$isoformSwitchAnalysis <- as.data.frame(cuffSplicing)
  }
  
  if( addIFmatrix ) {
    ifMatrix$isoform_id <- rownames(ifMatrix)
    rownames(ifMatrix) <- NULL
    ifMatrix <- ifMatrix[,c(
      which(colnames(ifMatrix) == 'isoform_id'),
      which(colnames(ifMatrix) != 'isoform_id')
    )]
    
    switchAnalyzeRlist$isoformRepIF <- ifMatrix
  }
  
  if (!quiet) {
    message("Done")
  }
  
  return(switchAnalyzeRlist)
}


# annotate AS for isofrom pair --------------------------------------------

extractSplicingEnrichmentNew <- function (pairwiseIsoComparison,localAS, splicingToAnalyze = "all",countGenes = FALSE, 
                                          returnResult = TRUE, plot = TRUE, minEventsForPlotting = 10,
                                          localTheme = theme_bw(base_size = 14),alpha = 0.05 ) {
  if (TRUE) {
    if (is.null(localAS)) {
      stop("The analsis of alternative must be performed before it can be summarized. Please use analyzeAlternativeSplicing() and try again.")
    }
    
    acceptedTypes <- c("A3", "A5", "ATSS", "ATTS", "ES", 
                       "IR", "MEE", "MES")
    if (!all(splicingToAnalyze %in% c("all", acceptedTypes))) {
      stop("The argument(s) supplied to 'typeOfconsequence' are not accepted. Please see ?analyzeAlternativeSplicing under details for description of which strings are allowed.")
    }
    splicingAnalyzed <- intersect(acceptedTypes, colnames(localAS))
    if ("all" %in% splicingToAnalyze) {
      splicingToAnalyze <- splicingAnalyzed
    }
    splicingNotAnalyzed <- setdiff(splicingToAnalyze, splicingAnalyzed)
    if (length(splicingNotAnalyzed)) {
      warning(paste("The following consequences appear not to have been analyzed and will therefor not be summarized:", 
                    paste(splicingNotAnalyzed, collapse = ", "), 
                    sep = " "))
    }
  }
  
  if (TRUE) {
    # localAS <- switchAnalyzeRlist$AlternativeSplicingAnalysis
    localAS <- localAS[which(localAS$isoform_id %in% pairwiseIsoComparison$isoformUpregulated | 
                               localAS$isoform_id %in% pairwiseIsoComparison$isoformDownregulated), 
                       ]
    m1 <- reshape2::melt(localAS[, c("isoform_id", "ES_genomic_start", 
                                     "MEE_genomic_start", "MES_genomic_start", "IR_genomic_start", 
                                     "A5_genomic_start", "A3_genomic_start", "ATSS_genomic_start", 
                                     "ATTS_genomic_start")], id.vars = "isoform_id")
    colnames(m1)[3] <- "genomic_start"
    m1$AStype <- sapply(strsplit(as.character(m1$variable), 
                                 "_"), function(x) x[1])
    m2 <- reshape2::melt(localAS[, c("isoform_id", "ES_genomic_end", 
                                     "MEE_genomic_end", "MES_genomic_end", "IR_genomic_end", 
                                     "A5_genomic_end", "A3_genomic_end", "ATSS_genomic_end", 
                                     "ATTS_genomic_end")], id.vars = "isoform_id")
    colnames(m2)[3] <- "genomic_end"
    m2$AStype <- sapply(strsplit(as.character(m2$variable), 
                                 "_"), function(x) x[1])
    localAS <- dplyr::inner_join(m1[, c("isoform_id", "AStype", 
                                        "genomic_start")], m2[, c("isoform_id", "AStype", 
                                                                  "genomic_end")], by = c("isoform_id", "AStype"))
    # localNMD <- data.frame(isoform_id = switchAnalyzeRlist$orfAnalysis$isoform_id, 
    #                        AStype = "NMD", genomic_start = ifelse(switchAnalyzeRlist$orfAnalysis$PTC, 
    #                                                               "0,0", "0"), genomic_end = ifelse(switchAnalyzeRlist$orfAnalysis$PTC, 
    #                                                                                                 "0,0", "0"), stringsAsFactors = FALSE)
    # localAS <- rbind(localAS, localNMD)
  }
  if (TRUE) {
    localConseq2 <- merge(pairwiseIsoComparison, localAS, 
                          by.x = "isoformUpregulated", by.y = "isoform_id")
    localConseq3 <- merge(localConseq2, localAS, by.x = c("isoformDownregulated", 
                                                          "AStype"), by.y = c("isoform_id", "AStype"), suffixes = c("_up", 
                                                                                                                    "_down"))
    localConseq3$genomic_start_up[which(is.na(localConseq3$genomic_start_up))] <- 0
    localConseq3$genomic_end_up[which(is.na(localConseq3$genomic_end_up))] <- 0
    localConseq3$genomic_start_down[which(is.na(localConseq3$genomic_start_down))] <- 0
    localConseq3$genomic_end_down[which(is.na(localConseq3$genomic_end_down))] <- 0
    localConseq3$coordinatsDifferent <- localConseq3$genomic_start_up != 
      localConseq3$genomic_start_down | localConseq3$genomic_end_up != 
      localConseq3$genomic_end_down
    localConseq4 <- localConseq3[which(localConseq3$coordinatsDifferent), 
                                 ]
    if (nrow(localConseq4) == 0) {
      stop("No alternative splicing differences were found")
    }
  }
  if (TRUE) {
    localConseq4 <- localConseq4[which(localConseq4$AStype %in% 
                                         splicingToAnalyze), ]
    if (!nrow(localConseq4)) {
      stop("No swithces with consequences were found")
    }
  }
  if (TRUE) {
    genomic_start_up <- strsplit(x = localConseq4$genomic_start_up, 
                                 split = ";")
    genomic_start_down <- strsplit(x = localConseq4$genomic_start_down, 
                                   split = ";")
    genomic_end_up <- strsplit(x = localConseq4$genomic_end_up, 
                               split = ";")
    genomic_end_down <- strsplit(x = localConseq4$genomic_end_down, 
                                 split = ";")
    localConseq4$nrGain <- 0
    localConseq4$nrLoss <- 0
    for (i in seq_along(genomic_start_up)) {
      localUp <- stringr::str_c(genomic_start_up[[i]], 
                                "_", genomic_end_up[[i]])
      localDn <- stringr::str_c(genomic_start_down[[i]], 
                                "_", genomic_end_down[[i]])
      localUp <- localUp[which(localUp != "0_0")]
      localDn <- localDn[which(localDn != "0_0")]
      localConseq4$nrGain[i] <- sum(!localUp %in% localDn)
      localConseq4$nrLoss[i] <- sum(!localDn %in% localUp)
    }
    terminiIndex <- which(localConseq4$AStype %in% c("ATSS", 
                                                     "ATTS"))
    localConseq4$nrGain[terminiIndex] <- 1 * sign(localConseq4$nrGain[terminiIndex])
    localConseq4$nrLoss[terminiIndex] <- 1 * sign(localConseq4$nrLoss[terminiIndex])
    localConseq4$nrDiff <- localConseq4$nrGain - localConseq4$nrLoss
  }
  #c("condition_1",   "condition_2", "AStype")
  gainLossBalance <- plyr::ddply(.data = localConseq4, .variables = c( "AStype"), .fun = function(aDF) {
    if (countGenes) {
      aDF2 <- plyr::ddply(aDF[, c("gene_id", "nrDiff")], 
                          .variables = "gene_id", function(aGene) {
                            data.frame(nrDiff = sum(aGene$nrDiff))
                          })
      localRes <- data.frame(nUp = sum(aDF2$nrDiff > 0), 
                             nDown = sum(aDF2$nrDiff < 0))
    }
    else {
      localRes <- data.frame(nUp = sum(aDF$nrDiff > 0), 
                             nDown = sum(aDF$nrDiff < 0))
    }
    if (localRes$nUp > 0 | localRes$nDown > 0) {
      localTest <- suppressWarnings(prop.test(localRes$nUp, 
                                              localRes$nUp + localRes$nDown))
      localRes$propUp <- localTest$estimate
      localRes$propUpCiLo <- min(localTest$conf.int)
      localRes$propUpCiHi <- max(localTest$conf.int)
      localRes$propUpPval <- localTest$p.value
    }
    else {
      localRes$propUp <- NA
      localRes$propUpCiLo <- NA
      localRes$propUpCiHi <- NA
      localRes$propUpPval <- NA
    }
    return(localRes)
  })
  
  
  
  if (TRUE) {
    gainLossBalance <- gainLossBalance[which(!is.na(gainLossBalance$propUp)), 
                                       ]
    gainLossBalance$propUpQval <- p.adjust(gainLossBalance$propUpPval, 
                                           method = "fdr")
    gainLossBalance$Significant <- gainLossBalance$propUpQval < 
      alpha
    gainLossBalance$Significant <- factor(gainLossBalance$Significant, 
                                          levels = c(FALSE, TRUE))
    # gainLossBalance$Comparison <- paste(gainLossBalance$condition_1, 
    #                                     "vs", gainLossBalance$condition_2, sep = "\n")
    gainLossBalance$Comparison <- paste("cond1","vs", "cond2", sep = "\n")
    myOrder <- plyr::ddply(gainLossBalance, .variables = "AStype", 
                           .fun = function(aDF) {
                             mean(aDF$propUp)
                           })
    myOrder <- myOrder$AStype[sort.list(myOrder$V1, decreasing = TRUE)]
    gainLossBalance$AStype <- factor(gainLossBalance$AStype, 
                                     levels = myOrder)
    gainLossBalance <- gainLossBalance[which((gainLossBalance$nUp + 
                                                gainLossBalance$nDown) >= minEventsForPlotting), 
                                       ]
    if (nrow(gainLossBalance) == 0) {
      stop("No features left for ploting after filtering with via \"minEventsForPlotting\" argument.")
    }
    gainLossBalance$nTot <- gainLossBalance$nUp + gainLossBalance$nDown
  }
  if (plot) {
    if (countGenes) {
      xText <- "Fraction of Switching Genes Primarly\nResulting in The Alternative Splicing Event Indicated\n(With 95% Confidence Interval)"
    } else {
      xText <- "Fraction of Switches Primarly\nResulting in Alternative Splicing Event Indicated\n(With 95% Confidence Interval)"
    }
    gainLossBalance$AStype2 <- paste(gainLossBalance$AStype, 
                                     "gain", "\n(paried with", gainLossBalance$AStype, 
                                     "loss)")
    g1 <- ggplot(data = gainLossBalance, aes(y = AStype2, 
                                             x = propUp, color = Significant)) + geom_point(aes(size = nTot)) + 
      geom_errorbarh(aes(xmax = propUpCiHi, xmin = propUpCiLo), 
                     height = 0.3) + facet_wrap(~Comparison) + geom_vline(xintercept = 0.5, 
                                                                          linetype = "dashed") + labs(x = xText, y = "Alternative Splicing Event\n(in isoform used more)") + 
      localTheme + theme(axis.text.x = element_text(angle = -45, 
                                                    hjust = 0, vjust = 1)) + scale_color_manual(name = paste0("FDR < ", 
                                                                                                              alpha), values = c("black", "red"), drop = FALSE) + 
      scale_size_continuous(name = "Observations") + guides(color = guide_legend(order = 1), 
                                                            size = guide_legend(order = 2)) + coord_cartesian(xlim = c(0, 
                                                                                                                       1))
    print(g1)
  }
  if (returnResult) {
    gainLossBalance$nTot <- NULL
    # localConseq5 <- localConseq3[, c("gene_id", "AStype", "isoformUpregulated", "isoformDownregulated", 
    #                                  "condition_1", "condition_2")]
    localConseq5 <- localConseq3[, c("gene_id", "AStype", "isoformUpregulated", "isoformDownregulated")]
    localConseq5$nrDiff <- localConseq4$nrDiff[match(stringr::str_c(localConseq5$isoformUpregulated, 
                                                                    localConseq5$isoformDownregulated, localConseq5$AStype), 
                                                     stringr::str_c(localConseq4$isoformUpregulated, localConseq4$isoformDownregulated, 
                                                                    localConseq4$AStype))]
    localConseq5$nrDiff[which(is.na(localConseq5$nrDiff))] <- 0
    localConseq5$iso_ref_up <- NULL
    localConseq5$iso_ref_down <- NULL
    localConseq5$ASchange <- dplyr::case_when(localConseq5$nrDiff > 
                                                0 ~ "Primarly gain", localConseq5$nrDiff < 0 ~ 
                                                "Primarly loss", TRUE ~ "No change")
    return(list(sumary=gainLossBalance,details=localConseq5))
    
  }
  
}


# import GTF --------------------------------------------------------------

importGTFNew <- function(
  pathToGTF,
  isoformNtFasta = NULL,
  extractAaSeq = FALSE,
  addAnnotatedORFs = TRUE,
  onlyConsiderFullORF = FALSE,
  removeNonConvensionalChr = FALSE,
  ignoreAfterBar = TRUE,
  ignoreAfterSpace = TRUE,
  ignoreAfterPeriod = FALSE,
  removeTECgenes = TRUE,
  PTCDistance = 50,
  quiet = FALSE
) {
  ### Test files
  if(TRUE) {
    ### Test existance of files
    if(TRUE) {
      if( pathToGTF == '' ) {
        stop(
          paste(
            'The \'pathToGTF\' argument does not lead anywhere (acutally you just suppled "" to the argument).',
            '\nDid you try to use the system.file("your/quant/dir/", package="IsoformSwitchAnalyzeR")',
            'to import your own data? The system.file() should only be used',
            'to access the example data stored in the IsoformSwitchAnalyzeR package.',
            'To access your own data simply provide the string to the directory with the data as:',
            '"path/to/quantification/".',
            sep=' '
          )
        )
      }
      if( ! file.exists(pathToGTF) ) {
        stop(
          paste(
            'The file pointed to with the \'pathToGTF\' argument does not exists.',
            '\nDid you accidentially make a spelling mistake or added a unwanted "/" infront of the text string?',
            sep=' '
          )
        )
      }
      
      if( !is.null(isoformNtFasta)) {
        if( !is.character( isoformNtFasta)) {
          stop('The \'isoformNtFasta\' argument must be a charachter string.')
        }
        
        if( any(isoformNtFasta == '') ) {
          stop(
            paste(
              'The \'isoformNtFasta\' argument does not lead anywhere (acutally you just suppled "" to the argument).',
              '\nDid you try to use the system.file("your/quant/dir/", package="IsoformSwitchAnalyzeR")',
              'to import your own data? The system.file() should only be used',
              'to access the example data stored in the IsoformSwitchAnalyzeR package.',
              'To access your own data simply provide the string to the directory with the data as:',
              '"path/to/quantification/".',
              sep=' '
            )
          )
        }
        
        if( any( ! sapply(isoformNtFasta, file.exists) ) ) {
          stop('At least one of the file(s) pointed to with \'isoformNtFasta\' seems not to exist.')
        }
        
        if( any(! grepl('\\.fa|\\.fasta|\\.fa.gz|\\.fasta.gz', isoformNtFasta)) ) {
          stop('The file pointed to via the \'isoformNtFasta\' argument does not seem to be a fasta file...')
        }
      }
      
      
    }
    
    # if( ! grepl('\\.gtf$|\\.gtf\\.gz$', pathToGTF, ignore.case = TRUE) ) {
    #   warning('The file pointed to by the "pathToGTF" argument appearts not to be a GTF file as it does not end with \'.gtf\' or \'.gtf.gz\' - are you sure it is the rigth file?')
    # }
    # 
  }
  
  # Read in from GTF file and create Rdata file for easy loading
  if (!quiet) {
    message('importing GTF (this may take a while)')
  }
  mfGTF <- rtracklayer::import(pathToGTF, format='gtf')
  
  ### Check GTF
  if (!all(c('transcript_id', 'gene_id') %in% colnames(mfGTF@elementMetadata))) {
    collumnsMissing <- paste(
      c('transcript_id', 'gene_id')[which(
        !c('transcript_id', 'gene_id') %in%
          colnames(mfGTF@elementMetadata)
      )], collapse = ', ')
    stop(
      paste(
        'The GTF file must contain the folliwing collumns',
        '\'transcript_id\', \'gene_id\' and.',
        collumnsMissing,
        'is missing.',
        sep = ' '
      )
    )
  }
  
  ### Reduce if nessesary
  if (removeNonConvensionalChr) {
    mfGTF <- mfGTF[which( ! grepl('_'  , as.character(mfGTF@seqnames))), ]
    mfGTF <- mfGTF[which( ! grepl('\\.', as.character(mfGTF@seqnames))), ]
    
    if (length(mfGTF) == 0) {
      stop('No exons were left after filtering',
           'with \'removeNonConvensionalChr\'.')
    }
    
    seqlevels(mfGTF) <- as.character(mfGTF@seqnames@values)
  }
  
  if( removeTECgenes ) {
    ### Ensembl
    if( 'gene_biotype' %in% colnames(mcols(mfGTF)) ) {
      mfGTF <- mfGTF[which(mfGTF$gene_biotype != 'TEC'),]
    }
    
    ### Gencode
    if( 'gene_type' %in% colnames(mcols(mfGTF)) ) {
      mfGTF <- mfGTF[which(mfGTF$gene_type != 'TEC'),]
    }
  }
  
  ### Potentially add version numbering
  if( TRUE ) {
    if(any( colnames(mcols(mfGTF)) == 'gene_version' )) {
      mfGTF$gene_id <- stringr::str_c(
        mfGTF$gene_id,
        '.',
        mfGTF$gene_version
      )
    }
    if(any( colnames(mcols(mfGTF)) == 'transcript_version' )) {
      mfGTF$transcript_id <- stringr::str_c(
        mfGTF$transcript_id,
        '.',
        mfGTF$transcript_version
      )
    }
  }
  
  ### Fix names
  if( ignoreAfterBar | ignoreAfterSpace | ignoreAfterPeriod) {
    
    mfGTF$transcript_id <- IsoformSwitchAnalyzeR:::fixNames(
      nameVec = mfGTF$transcript_id,
      ignoreAfterBar = ignoreAfterBar,
      ignoreAfterSpace = ignoreAfterSpace,
      ignoreAfterPeriod = ignoreAfterPeriod
    )
    
  }
  
  ### Make annoation
  if(TRUE) {
    if (!quiet) {
      message('converting GTF to switchAnalyzeRlist')
    }
    exonAnoationIndex <- which(mfGTF$type == 'exon')
    
    colsToExtract <- c(
      'transcript_id', 'gene_id', 'gene_name',
      'gene_type','gene_biotype',     # respectively gencode and ensembl gene type col
      'transcript_biotype','transcript_type'
    )
    myIso <-
      as.data.frame(unique(mfGTF@elementMetadata[
        exonAnoationIndex,
        na.omit(match(colsToExtract, colnames(mfGTF@elementMetadata)))]
      ))
    
    ### Handle columns not extracted
    if (is.null(myIso$gene_name)) {
      myIso$gene_name <- NA
    }
    
    ### Handle columns with multiple options
    geneTypeCol <- which(colnames(myIso) %in% c('gene_type','gene_biotype'))
    if( length(geneTypeCol) == 0 ) {
      myIso$geneType <- NA
    } else {
      myIso$geneType <- myIso[,geneTypeCol]
    }
    
    isoTypeCol <- which(colnames(myIso) %in% c('transcript_biotype','transcript_type'))
    if( length(geneTypeCol) == 0 ) {
      myIso$isoType <- NA
    } else {
      myIso$isoType <- myIso[,isoTypeCol]
    }
    
    ### Make annotation
    myIsoAnot <- data.frame(
      isoform_id = myIso$transcript_id,
      gene_id = myIso$gene_id,
      condition_1 = "plaseholder1",
      condition_2 = "plaseholder2",
      gene_name = myIso$gene_name,
      gene_biotype = myIso$geneType,
      iso_biotype = myIso$isoType,
      class_code = '=',
      gene_overall_mean = 0,
      gene_value_1 = 0,
      gene_value_2 = 0,
      gene_stderr_1 = NA,
      gene_stderr_2 = NA,
      gene_log2_fold_change = 0,
      gene_p_value = 1,
      gene_q_value = 1,
      iso_overall_mean = 0,
      iso_value_1 = 0,
      iso_value_2 = 0,
      iso_stderr_1 = NA,
      iso_stderr_2 = NA,
      iso_log2_fold_change = 0,
      iso_p_value = 1,
      iso_q_value = 1,
      IF_overall = NA,
      IF1 = NA,
      IF2 = NA,
      dIF = NA,
      isoform_switch_q_value = NA,
      gene_switch_q_value = NA,
      stringsAsFactors = FALSE
    )
    
  }
  
  ### Add CDS annoation from GTF file inc convertion to transcript coordinats
  if (addAnnotatedORFs) {
    # test whether any CDS are found
    if (any(mfGTF$type == 'CDS')) {
      if (!quiet) {
        message('converting annotated CDSs')
      }
      myCDS <-
        sort(mfGTF[which(mfGTF$type == 'CDS'), 'transcript_id'])
      myCDSedges <-
        suppressMessages(unlist(range(
          split(myCDS[, 0], f = myCDS$transcript_id)
        )))  # Extract EDGEs
      myCDSedges$id <- names(myCDSedges)
      names(myCDSedges) <- NULL
      
      if (onlyConsiderFullORF) {
        fullyAnnoated <-
          as.data.frame(sort(
            mfGTF[which(
              mfGTF$type %in% c('start_codon', 'stop_codon')),
              c('transcript_id', 'type')]))
        fullyAnnoatedSplit <-
          split(as.character(fullyAnnoated$type),
                f = fullyAnnoated$transcript_id)
        fullyAnnoatedCount <-
          sapply(fullyAnnoatedSplit, function(x)
            length(unique(x)))
        toKeep <-
          names(fullyAnnoatedCount[which(fullyAnnoatedCount == 2)])
        
        
        myCDSedges <-
          myCDSedges[which(myCDSedges$id %in% toKeep), ]
      }
      
      ### Extract Exons
      localExons <- mfGTF[exonAnoationIndex, 'transcript_id']
      localExons <-
        localExons[which(
          as.character(localExons@strand) %in% c('+', '-')), ]
      localExons <-
        localExons[which(localExons$transcript_id %in% myCDSedges$id), ]
      
      localExons <-
        localExons[order(localExons$transcript_id,
                         start(localExons),
                         end(localExons)), ]
      localExons$exon_id <-
        paste('exon_', 1:length(localExons), sep = '')
      
      ### Extract strand specific ORF info
      cds <- as.data.frame(myCDSedges)
      # start
      plusIndex <- which(cds$strand == '+')
      annoatedStartGRangesPlus <-
        GRanges(
          cds$seqnames[plusIndex],
          IRanges(
            start = cds$start[plusIndex],
            end = cds$start[plusIndex]),
          strand = cds$strand[plusIndex],
          id = cds$id[plusIndex]
        )
      minusIndex <- which(cds$strand == '-')
      annoatedStartGRangesMinus <-
        GRanges(
          cds$seqnames[minusIndex],
          IRanges(
            start = cds$end[minusIndex],
            end = cds$end[minusIndex]),
          strand = cds$strand[minusIndex],
          id = cds$id[minusIndex]
        )
      
      annoatedStartGRanges <-
        c(annoatedStartGRangesPlus,
          annoatedStartGRangesMinus)
      annoatedStartGRanges$orf_id <-
        paste('cds_', 1:length(annoatedStartGRanges), sep = '')
      
      # end
      annoatedEndGRangesPlus  <-
        GRanges(
          cds$seqnames[plusIndex],
          IRanges(
            start = cds$end[plusIndex],
            end = cds$end[plusIndex]),
          strand = cds$strand[plusIndex],
          id = cds$id[plusIndex]
        )
      annoatedEndGRangesMinus <-
        GRanges(
          cds$seqnames[minusIndex],
          IRanges(
            start = cds$start[minusIndex],
            end = cds$start[minusIndex]),
          strand = cds$strand[minusIndex],
          id = cds$id[minusIndex]
        )
      
      annoatedEndGRanges <-
        c(annoatedEndGRangesPlus, annoatedEndGRangesMinus)
      annoatedEndGRanges$orf_id <-
        paste('stop_', 1:length(annoatedEndGRanges), sep = '')
      
      # combine
      annotatedORFGR <-
        c(annoatedStartGRanges, annoatedEndGRanges)
      
      
      ### Idenetify overlapping CDS and exons as well as the annoate transcript id
      suppressWarnings(overlappingAnnotStart <-
                         as.data.frame(
                           findOverlaps(
                             query = localExons,
                             subject = annotatedORFGR,
                             ignore.strand = FALSE
                           )
                         ))
      if (!nrow(overlappingAnnotStart)) {
        stop(
          'No overlap between CDS and transcripts were found. This is most likely due to a annoation problem around chromosome name.'
        )
      }
      
      # Annoate overlap ids
      overlappingAnnotStart$transcript_id <-
        localExons$transcript_id[overlappingAnnotStart$queryHits]
      overlappingAnnotStart$exon_id <- localExons$exon_id[
        overlappingAnnotStart$queryHits
        ]
      
      overlappingAnnotStart$cdsTranscriptID <- annotatedORFGR$id[
        overlappingAnnotStart$subjectHits
        ]
      overlappingAnnotStart$orf_id <- annotatedORFGR$orf_id[
        overlappingAnnotStart$subjectHits
        ]
      
      # subset to annoateted overlap
      overlappingAnnotStart <-
        overlappingAnnotStart[which(
          overlappingAnnotStart$transcript_id ==
            overlappingAnnotStart$cdsTranscriptID
        ), c('transcript_id',
             'exon_id',
             'cdsTranscriptID',
             'orf_id')]
      
      # annoate with genomic site
      overlappingAnnotStart$orfGenomic <-
        start(annotatedORFGR)[match(
          overlappingAnnotStart$orf_id, annotatedORFGR$orf_id
        )]
      
      
      ### Enrich exon information
      myExons <-
        as.data.frame(localExons[which(
          localExons$transcript_id %in%
            overlappingAnnotStart$transcript_id),])
      
      # Strand
      myExonPlus <- myExons[which(myExons$strand == '+'), ]
      myExonMinus <- myExons[which(myExons$strand == '-'), ]
      
      plusSplit <-
        split(myExonPlus$width, myExonPlus$transcript_id)
      minusSplit <-
        split(myExonMinus$width, myExonMinus$transcript_id)
      
      # cumsum
      myExonPlus$cumSum <-
        unlist(sapply(plusSplit , function(aVec) {
          cumsum(c(0, aVec))[1:(length(aVec))]
        }))
      myExonMinus$cumSum <-
        unlist(sapply(minusSplit, function(aVec) {
          cumsum(c(0, rev(aVec)))[(length(aVec)):1] # reverse
        }))
      
      # exon number
      myExonPlus$nrExon <-
        unlist(sapply(plusSplit, function(aVec) {
          1:length(aVec)
        }))
      myExonMinus$nrExon <-
        unlist(sapply(minusSplit, function(aVec) {
          1:length(aVec)
        }))
      
      # total nr exons
      myExonPlus$lastExonIndex <-
        unlist(sapply(plusSplit, function(aVec) {
          rep(length(aVec), length(aVec))
        }))
      myExonMinus$lastExonIndex <-
        unlist(sapply(minusSplit, function(aVec) {
          rep(1, length(aVec))
        }))
      
      # final exon exon junction trancipt position
      myExonPlus$finalJunctionPos <-
        unlist(sapply(plusSplit , function(aVec) {
          rep(cumsum(c(0, aVec))[length(aVec)], times = length(aVec))
        }))
      myExonMinus$finalJunctionPos <-
        unlist(sapply(minusSplit, function(aVec) {
          rep(cumsum(c(0, rev(
            aVec
          )))[length(aVec)], times = length(aVec))
        }))
      
      myExons2 <- rbind(myExonPlus, myExonMinus)
      
      ### Annoate with exon information
      matchIndex <-
        match(overlappingAnnotStart$exon_id, myExons2$exon_id)
      overlappingAnnotStart$strand <- myExons2$strand[matchIndex]
      overlappingAnnotStart$exon_start <- myExons2$start[matchIndex]
      overlappingAnnotStart$exon_end <- myExons2$end[matchIndex]
      overlappingAnnotStart$exon_cumsum <- myExons2$cumSum[matchIndex]
      overlappingAnnotStart$exon_nr <- myExons2$nrExon[matchIndex]
      overlappingAnnotStart$lastExonIndex <-
        myExons2$lastExonIndex[matchIndex]
      overlappingAnnotStart$finalJunctionPos <-
        myExons2$finalJunctionPos[matchIndex]
      
      ### Annoate with transcript coordinats
      overlappingAnnotStartPlus <-
        overlappingAnnotStart[which(
          overlappingAnnotStart$strand == '+'), ]
      overlappingAnnotStartPlus$orfTranscript <-
        overlappingAnnotStartPlus$exon_cumsum + (
          overlappingAnnotStartPlus$orfGenomic -
            overlappingAnnotStartPlus$exon_start
        ) + 1
      overlappingAnnotStartPlus$junctionDistance <-
        overlappingAnnotStartPlus$finalJunctionPos -
        overlappingAnnotStartPlus$orfTranscript + 3 # +3 because the ORF does not include the stop codon - but it should in this calculation
      
      overlappingAnnotStartMinus <-
        overlappingAnnotStart[which(
          overlappingAnnotStart$strand == '-'), ]
      overlappingAnnotStartMinus$orfTranscript <-
        overlappingAnnotStartMinus$exon_cumsum + (
          overlappingAnnotStartMinus$exon_end -
            overlappingAnnotStartMinus$orfGenomic
        ) + 1
      overlappingAnnotStartMinus$junctionDistance <-
        overlappingAnnotStartMinus$finalJunctionPos -
        overlappingAnnotStartMinus$orfTranscript + 3 # +3 because the ORF does not include the stop codon - but it should in this calculation
      
      overlappingAnnotStart2 <-
        rbind(overlappingAnnotStartPlus,
              overlappingAnnotStartMinus)
      overlappingAnnotStart2 <-
        overlappingAnnotStart2[order(
          overlappingAnnotStart2$transcript_id,
          overlappingAnnotStart2$exon_start,
          overlappingAnnotStart2$exon_end
        ), ]
      
      ### devide into start and stop
      starInfo <-
        overlappingAnnotStart2[which(
          grepl('^cds', overlappingAnnotStart2$orf_id)), ]
      stopInfo <-
        overlappingAnnotStart2[which(
          grepl('^stop', overlappingAnnotStart2$orf_id)), ]
      
      ### predict PTC
      stopInfo$PTC <-
        stopInfo$exon_nr != stopInfo$lastExonIndex &
        stopInfo$junctionDistance > PTCDistance
      
      ### Merge the data
      starInfo2 <-
        unique(starInfo[, c('transcript_id',
                            'orfGenomic',
                            'exon_nr',
                            'orfTranscript')])
      colnames(starInfo2) <-
        c('isoform_id',
          'orfStartGenomic',
          'orfStarExon',
          'orfTransciptStart')
      
      stopInfo2 <-
        unique(stopInfo[, c(
          'transcript_id',
          'orfGenomic',
          'exon_nr',
          'orfTranscript',
          'junctionDistance',
          'lastExonIndex',
          'PTC'
        )])
      colnames(stopInfo2) <-
        c(
          'isoform_id',
          'orfEndGenomic',
          'orfEndExon',
          'orfTransciptEnd',
          'stopDistanceToLastJunction',
          'stopIndex',
          'PTC'
        )
      
      orfInfo <- dplyr::inner_join(starInfo2, stopInfo2, by = 'isoform_id')
      orfInfo$orfTransciptLength  <-
        orfInfo$orfTransciptEnd - orfInfo$orfTransciptStart + 1
      
      # reorder
      orfInfo <-
        orfInfo[, c(
          'isoform_id',
          'orfTransciptStart',
          'orfTransciptEnd',
          'orfTransciptLength',
          'orfStarExon',
          'orfEndExon',
          'orfStartGenomic',
          'orfEndGenomic',
          'stopDistanceToLastJunction',
          'stopIndex',
          'PTC'
        )]
      
      # make sure all ORFs are annotated (with NAs)
      orfInfo <-
        dplyr::full_join(orfInfo,
                         unique(myIsoAnot[, 'isoform_id', drop = FALSE]),
                         by = 'isoform_id',
                         all = TRUE)
      
    } else {
      # if no CDS were found
      warning(paste(
        'No CDS was found in the GTF file. Please make sure the GTF',
        'file have the CDS "feature" annotation. Adding NAs instead'
      ))
      
      orfInfo <- data.frame(
        isoform_id = unique(myIsoAnot$isoform_id),
        orfTransciptStart = NA,
        orfTransciptEnd = NA,
        orfTransciptLength = NA,
        orfStarExon = NA,
        orfEndExon = NA,
        orfStartGenomic = NA,
        orfEndGenomic = NA,
        stopDistanceToLastJunction = NA,
        stopIndex = NA,
        PTC = NA,
        stringsAsFactors = FALSE
      )
    }
    
    ### add to iso annotation
    myIsoAnot$PTC <-
      orfInfo$PTC[match(myIsoAnot$isoform_id, orfInfo$isoform_id)]
    
  }
  
  ### Handle sequence input
  if(TRUE) {
    addIsoformNt <- FALSE
    
    if( !is.null(isoformNtFasta) ) {
      isoformNtSeq <- do.call(
        c,
        lapply(isoformNtFasta, function(aFile) {
          Biostrings::readDNAStringSet(
            filepath = isoformNtFasta, format = 'fasta'
          )
        })
      )
      
      if(!is(isoformNtSeq, "DNAStringSet")) {
        stop('The fasta file supplied to \'isoformNtFasta\' does not contain the nucleotide (DNA) sequence...')
      }
      
      ### Fix names
      if( ignoreAfterBar | ignoreAfterSpace | ignoreAfterPeriod) {
        
        names(isoformNtSeq) <- IsoformSwitchAnalyzeR:::fixNames(
          nameVec = names(isoformNtSeq),
          ignoreAfterBar = ignoreAfterBar,
          ignoreAfterSpace = ignoreAfterSpace,
          ignoreAfterPeriod = ignoreAfterPeriod
        )
      }
      
      ### Subset to annotated isoforms
      isoformNtSeq <- isoformNtSeq[which(
        names(isoformNtSeq) %in% myIsoAnot$isoform_id
      )]
      
      ### Remove potential name duplication
      isoformNtSeq <- isoformNtSeq[which(
        ! duplicated(names(isoformNtSeq))
      )]
      
      if( ! all(myIsoAnot$isoform_id %in% names(isoformNtSeq)) ) {
        warning(
          paste(
            'The fasta file supplied to \'isoformNtFasta\' does not contain the',
            'nucleotide (DNA) sequence for all isoforms annotated and will not be added',
            '\nSpecifically:\n',
            length(unique(myIsoAnot$isoform_id)), 'isoforms were annotated in the GTF\n',
            length(unique(names(isoformNtSeq))), 'isoforms have a sequence.\n',
            'Only', length(intersect(names(isoformNtSeq), myIsoAnot$isoform_id)), 'overlap.\n',
            length(setdiff(unique(myIsoAnot$isoform_id), names(isoformNtSeq))), 'annoated isoforms isoforms had no corresponding nucleotide sequence\n',
            
            '\nIf there is no overlap (as in zero or close) there are two options:\n',
            '1) The files do not fit together (different databases, versions etc)',
            '(no fix except using propperly paired files).\n',
            '2) It is somthing to do with how the isoform ids are stored in the different files.',
            'This problem might be solvable using some of the',
            '\'ignoreAfterBar\', \'ignoreAfterSpace\' or \'ignoreAfterPeriod\' arguments.\n',
            '    3 Examples from GTF are :',
            paste0( sample(unique(myIsoAnot$isoform_id), 3), collapse = ', '),'\n',
            '    3 Examples of isoform sequence are  :',
            paste0( sample(names(isoformNtSeq), 3), collapse = ', '),'\n',
            
            
            '\nIf there is a large overlap but still far from complete there are 3 possibilites:\n',
            '1) The files do not fit together (different databases versions)',
            '(no fix except using propperly paired files).\n',
            '2) The isoforms quantified have their nucleotide sequence stored in multiple fasta files (common for Ensembl).',
            'Just supply a vector with the path to each of them to the \'isoformNtFasta\' argument.\n',
            '3) One file could contain non-chanonical chromosomes while the other do not',
            '(might be solved using the \'removeNonConvensionalChr\' argument.)\n\n',
            sep = ' '
          )
        )
      } else {
        addIsoformNt <- TRUE
      }
    }
  }
  
  # Create exon_features grange
  myExons <-
    sort(mfGTF[exonAnoationIndex , c('transcript_id', 'gene_id')])
  colnames(myExons@elementMetadata) <- c('isoform_id', 'gene_id')
  
  # Collaps ajecent exons (without any intron between)
  if(TRUE) {
    ### Reduce ajecent exons
    tmp <- unlist(
      reduce(
        split(
          myExons,
          myExons$isoform_id
        )
      )
    )
    ### Add isoform id
    tmp$isoform_id <- tmp@ranges@NAMES
    tmp@ranges@NAMES <- NULL
    
    ### add gene id
    tmp$gene_id <-myExons$gene_id[match(
      tmp$isoform_id, myExons$isoform_id
    )]
    
    ### sort
    tmp <- tmp[sort.list(tmp$isoform_id),]
    
    ### Overwrite
    myExons <- tmp
  }
  
  # create replicates
  nrRep <-
    data.frame(
      condition = c('plaseholder1', 'plaseholder2'),
      nrReplicates = c(NA, NA),
      row.names = NULL,
      stringsAsFactors = FALSE
    )
  
  # create dummy feature
  repExp <- data.frame(
    isoform_id = myIsoAnot$isoform_id,
    plaseholder1 = NA,
    plaseholder2 = NA,
    stringsAsFactors = FALSE
  )
  
  designMatrix <-
    data.frame(
      sampleID = c('plaseholder1', 'plaseholder2'),
      condition = c('plaseholder1', 'plaseholder2'),
      stringsAsFactors = FALSE
    )
  
  ### Create switchList
  localSwichList <- createSwitchAnalyzeRlist(
    isoformFeatures = myIsoAnot,
    exons = myExons,
    designMatrix = designMatrix,
    isoformCountMatrix = repExp,
    sourceId = 'gtf'
  )
  
  if (addAnnotatedORFs) {
    # subset to those in list
    orfInfo <-
      orfInfo[which(orfInfo$isoform_id %in%
                      localSwichList$isoformFeatures$isoform_id), ]
    
    # check for negative ORF lengths
    isoformsToRemove <-
      orfInfo$isoform_id[which(orfInfo$orfTransciptLength < 0)]
    if (length(isoformsToRemove)) {
      genesToRemove <-
        localSwichList$isoformFeatures$gene_id[which(
          localSwichList$isoformFeatures$isoform_id %in%
            isoformsToRemove)]
      localSwichList <-
        subsetSwitchAnalyzeRlist(
          localSwichList,
          !localSwichList$isoformFeatures$gene_id %in% genesToRemove
        )
      
      warning(
        paste(
          length(genesToRemove),
          'genes where removed due to negative ORF lengths. This',
          'typically occures because gene_id are not unique',
          '(meaning are found multiple places accorss the genome).',
          'Please note there might still be duplicated gene_id',
          'located on the same chromosome.',
          sep = ' '
        )
      )
    }
    localSwichList$orfAnalysis <- orfInfo[which(
      orfInfo$isoform_id %in% localSwichList$isoformFeatures$isoform_id)
      ,]
  }
  
  ### Add nucleotide sequence
  if(addIsoformNt) {
    localSwichList$ntSequence <- isoformNtSeq[which(
      names(isoformNtSeq) %in% localSwichList$isoformFeatures$isoform_id
    )]
    
    if(addAnnotatedORFs & extractAaSeq) {
      localSwichList <- extractSequence(
        switchAnalyzeRlist = localSwichList,
        onlySwitchingGenes = FALSE,
        writeToFile = FALSE,
        extractNTseq = TRUE,
        extractAAseq = TRUE,
        addToSwitchAnalyzeRlist = TRUE
      )
    }
  }
  
  
  # Return SpliceRList
  return(localSwichList)
}
