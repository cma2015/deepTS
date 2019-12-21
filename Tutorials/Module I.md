**Module I: Generation of transcript expression matrix**

<table>
    <tr>
      <td font-weight:bold>Functions</td>
      <td font-weight:bold>Applications</td>
      <td font-weight:bold>Input files</td>
      <td font-weight:bold>Main output files</td>
      <td font-weight:bold>Programs</td>
     </tr>
     <tr>
      <td>Preparation of high-quality RNA-Seq data</td>
      <td>Upload or retrieve raw RNA-Seq data and trim low-quality reads and sequencing adapters</td>
      <td>Raw RNA-Seq data or accession IDs/ftp address of data</td>
      <td>High-quality RNA-Seq data</td>
      <td>[fastp](https://academic.oup.com/bioinformatics/article/34/17/i884/5093234) (version 0.20.0)</td>
     </tr>
     <tr>
      <td rowspan="9">Construction of the transcriptome map</td>
      <td rowspan="9">Construct transcriptome map by combining reference, assembled and/or PacBio transcriptome</td>
      <td rowspan="9">Reference genome sequence; reference genome annotation; high-quality RNA-Seq data (FLNC PacBio reads)</td>
      <td rowspan="9">Transcriptome map; mapping results (BAM files)</td>
      <td>HISAT (version 2.1.0)(<https://www.nature.com/articles/nmeth.3317>)</td>
     </tr>
     <tr>
      <td>SAMTools (version 1.10)(<https://academic.oup.com/bioinformatics/article/25/16/2078/204688>)</td>
     </tr>
     <tr>
      <td>BEDTools (version 2.29.0)(<https://academic.oup.com/bioinformatics/article/26/6/841/244688>)</td>
     </tr>
     <tr>     
      <td>StringTie (version 1.3.4)(<https://www.nature.com/articles/nbt.3122>)</td>
     </tr>
     <tr>  
      <td>Cufflinks (version 2.2.1)(<https://www.nature.com/articles/nbt.1621>)</td>
     </tr>
     <tr>  
      <td>CPC2 (version 0.1)(<https://academic.oup.com/nar/article/45/W1/W12/3831091>)</td>
     </tr>
     <tr>  
      <td>DIAMOND (version 0.9.29)(<https://www.nature.com/articles/nmeth.3176>)</td>
     </tr>
     <tr>  
      <td>featureCounts (version 2.0.0)(<https://academic.oup.com/bioinformatics/article/30/7/923/232889>)</td>
     </tr>
     <tr>  
      <td>GMAP (version  2015-09-29)(<https://academic.oup.com/bioinformatics/article/21/9/1859/409207>)</td>
     </tr>
     <tr>
      <td rowspan="2">Generation of expression matrix</td>
      <td rowspan="2">Estimate expression abundance of genes and transcripts in terms of TPM</td>
      <td rowspan="2">Transcriptome map; mapping results (BAM files)</td>
      <td rowspan="2">Expression matrix</td>
      <td>StringTie (version 1.3.4)(<https://www.nature.com/articles/nbt.3122>)</td>
     </tr>
     <tr>  
      <td>sva (version 3.34.0)(<https://academic.oup.com/bioinformatics/article/28/6/882/311263>)</td>
     </tr>
    </table>
