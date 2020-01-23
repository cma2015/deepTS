# **Module I: Generation of transcript expression matrix**

<table class="fl-table">
    <thead>
    <tr>
        <th width="15%">Functional modules</th>
        <th width="15%">Funtions</th>
        <th width="20%">Applications</th>
        <th width="20%">Input Files</th>
        <th width="15%">Main Output Files</th>
  <th width="15%">Programs</th>
    </tr>
    </thead>
    <tbody>
    <tr>
        <td rowspan="4">Identification of TS events</td>
    </tr>
    <tr>
        <td>Pairwise transcriptome comparison</td>
        <td>TS events identifiacton for pairwise transcriptome comparison</td>
        <td>Transcript expression abundance matrix; Gene-transcript ID matrix; Experiment design matrix</td>
        <td>TS events matrix</td>
  <td><a href="https://github.com/comprna/SUPPA">supppa2</a></td>
    </tr>
    <tr>
        <td>Time-series transcriptome comparison</td>
        <td>TS events identifiacton for time-series transcriptome comparison</td>
        <td>Transcript expression abundance matrix; Gene-transcript ID matrix; Experiment design matrix</td>
        <td>TS events matrix</td>	
  <td><a href="https://github.com/wyguo/TSIS">TSIS</a></td>
    </tr>
  <tr>
        <td>Population transcriptome comparison</td>
        <td>TS events identifiacton for population transcriptome comparison</td>
        <td>Genotype matrix, Transcript expression abundance matrix; Gene-transcript ID mapping matrix</td>
        <td>GWAS and TS events result</td>	
  <td><a href="http://www.zzlab.net/GAPIT/">GAPIT</a></td>
    </tr>
    <tbody>
</table>

# **Example data**

The example data for mudule II is can be download here(https://github.com/cma2015/DeepTS/blob/master/Test_data/Module2)

## Pairwise transcriptome comparison
- **Genome annotation file** is in the directory named `pairwiseTS_annotation_file.gtf`
- **Transcript expression abundance matrix** is in the directorynamed `pairwiseTS_transcript_expression.tsv`
- **Experiment design matrix** is in the directory named `pairwiseTS_experiment_design.tsv`

- **Gene-transcript ID matrix** is in the directory named `gene_transcriptID.tsv` 
The first, second and third columns are transcripts, genes and corresponding chromosomes. Chromosomes are in digital format.

## Time-series transcriptome comparison

- **Transcript expression abundance matrix** is in the directory named `timeSeriesTS_transcript_expression.tsv`
- **Experiment design matrix*** is in the directory named `timeSeriesTS_experiment_design.tsv`
- **Gene-transcript ID matrix** is in the directory named `gene_transcriptID.tsv`

## Population transcriptome comparison
- **Genotype data** is in the directory named `popultaionTS_genotype_hmpfile.tsv`
the chromsome names in hmp file are in digital format. 
- **Transcript expression abundance matrix** is in the directory named `popultaionTS_transcript_expression.tsv`
- **Gene-transcript ID matrix** is in the directory named `gene_transcriptID.tsv`

