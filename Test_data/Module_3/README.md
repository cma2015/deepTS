# [Exploration of TS events  ](https://github.com/cma2015/DeepTS/blob/master/Tutorials/Test_data/Module3)

## User-defined filtering
- **Inital TS events or significant SNP Matrix** is in the directory named `Inital_time-series_TS_events.tsv` These files come from the tools in modlue II **Identification of TS events events**.


## Characteristics analysis

- **TS events Matrix** is in the directory named `Filtered_time-series_TS_events.tsv`. 
- **Genome feature Matrix** is in the directory named `feature_genome_file.tsv`.
First to third columns are genomic coordinates, other columns are different characteristics.
- **Chromosome length** is in the directory` named `chromosome_length.tsv`
- **Gene feature flie** is in the directory named `feature_gene_file.tsv`
The first column are gene IDs, the other columns are different characteristics. Features with the suffix ".den" will be plotted in the density diagram. Features with the suffix ".pie" will be plotted in the piechart diagram.
- **Transcript feature flie** is in the directory named `feature_transcript_file.tsv`. The first column are transcript IDs, the other columns are different characteristics. Features with the suffix ".den" will be plotted in the density diagram. Features with the suffix ".pie" will be plotted in the piechart diagram.
- **Protein feature flie** is in the directory named `feature_protein_file.tsv`. The first column are transcript IDs, the other columns are different characteristics. Features with the suffix ".den" will be plotted in the density diagram. Features with the suffix ".pie" will be plotted in the piechart diagram.


## Expression pattern visualization

- **TS events Matrix** is in the directory named `Filtered_time-series_TS_events.tsv` The first column must be the group name representing the results of the different TS events. 

- **sqlitedb** is in the directory named `sqtlite_time-series_exploration_analysis.sqlite` created by **Creat SQLite file** in **other tools**.

## Venn diagram

- **TS events Matrix** is in the directory named `Filtered_time-series_TS_events.tsv` The first column must be the group name representing the results of the different TS events. 

## GO enrichment analysis

- **TS events Matrix** is in the directory named `Filtered_time-series_TS_events.tsv`.
- **GO terms ID** is in the directory named `goterm_test.tsv`.

## Alternative splicing

- **TS events Matrix** is in the directory named `Filtered_time-series_TS_events.tsv`.
- **Genome annotation file** is in the directory named `exploreTS_annotation_file.gtf`.
