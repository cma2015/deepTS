# **Module I: Generation of transcript expression matrix**

  <table class="fl-table">
  <thead>
    <tr>
      <th width="15%">Functional modules</th>
      <th width="15%">Functions</th>
      <th width="20%">Applications</th>
      <th width="20%">Input files</th>
      <th width="15%">Main output files</th>
      <th width="15%">Programs</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td rowspan="16">Generation of transcript
      expression matrix</td>
    </tr>
    <tr>
      <td>Preparation of high-quality
      RNA-Seq data</td>
      <td>Upload or retrieve raw RNA-Seq
      data and trim low-quality reads
      and sequencing adapters</td>
      <td>Raw RNA-Seq data or accession
      IDs/ftp address of data</td>
      <td>High-quality RNA-Seq data</td>
      <td><a href="https://academic.oup.com/bioinformatics/article/34/17/i884/5093234">fastp (version 0.20.0)</a></td>
    </tr>
    <tr>
      <td rowspan="12">Construction of the
      transcriptome map</td>
    </tr>
    <tr>
      <td rowspan="11">Construct transcriptome map
      by combining reference, assembled and/or 
      PacBio transcriptome</td>
    </tr>
    <tr>
      <td rowspan="10">Reference genome sequence; 
      reference genome annotation;high-quality 
      RNA-Seq data(FLNC PacBio reads)</td>
    </tr>
    <tr>
      <td rowspan="9">Transcriptome map;</br>
      mapping results (BAM files)</td>						
      <td><a href="https://www.nature.com/articles/nmeth.3317">HISAT (version 2.1.0)</a></td>					
    </tr>
    <tr>
      <td><a href="https://academic.oup.com/bioinformatics/article/25/16/2078/204688">SAMTools (version 1.10)</a></td>
    </tr>
    <tr>
      <td><a href="https://academic.oup.com/bioinformatics/article/26/6/841/244688">BEDTools (version 2.29.0)</td>	
    </tr>
    <tr>
      <td><a href="https://www.nature.com/articles/nbt.3122">StringTie (version 1.3.4)</td>				
    </tr>
    <tr>
      <td><a href="https://www.nature.com/articles/nbt.1621">Cufflinks (version 2.2.1)</td>
    </tr>
    <tr>
      <td><a href="https://academic.oup.com/nar/article/45/W1/W12/3831091">CPC2 (version 0.1)</td>
    </tr>
    <tr>
      <td><a href="https://www.nature.com/articles/nmeth.3176">DIAMOND (0.9.29)</td>						
    </tr>
    <tr>
      <td><a href="https://academic.oup.com/bioinformatics/article/30/7/923/232889">featureCounts (version 2.0.0)</td>				
    </tr>
    <tr>
      <td><a href="https://academic.oup.com/bioinformatics/article/21/9/1859/409207">GMAP (version  2015-09-29)</td>					
    </tr>			
    <tr>
        <td rowspan="2">Generation of expression matrix</td>
      <td rowspan="2">Estimate expression abundance
      of genes and transcripts 
      in terms of TPM</td>	
      <td rowspan="2">Transcriptome map;
      mapping results (BAM files)</td>	
      <td rowspan="2">Expression matrix</td>	
      <td><a href="https://www.nature.com/articles/nbt.3122">StringTie (version 1.3.4)</td>						
    </tr>
    <tr>		
      <td><a href="https://academic.oup.com/bioinformatics/article/28/6/882/311263">sva (version 3.34.0)</td>
    </tr>
  <tbody>
  </table>