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
            <td rowspan="8">Exploration of TS events</td>
        </tr>
        <tr>
            <td>User-defined filtering</td>
            <td>Filter three types of TS events</td>
            <td>Inital TS events matrix</td>
            <td>Filtered TS events matrix</td>
			<td>In-house scripts</td>
        </tr>
        <tr>
            <td>Characteristics analysis</td>
            <td>Explore the characteristics of TS events at the genome, gene, transcript, and protein levels</td>
            <td>TS events matrix and feature matrix</td>
            <td>html report</td>	
			<td>In-house scripts</td>
        </tr>
		 <tr>
            <td>Expression pattern visualization</td>
            <td>Display expression pattern, exon-intron structure and protein domain of transcripts in TS events</td>
            <td>sQlite file,TS events matrix and transcript ID </td>
            <td>PDF file</td>	
			<td>In-house scripts</td>
        </tr>
		 <tr>
            <td>Venn diagram</td>
            <td>Display similatities, differences and enrichment for different group TS events matrix</td>
            <td>TS events matrix</td>
            <td>PDF file</td>	
			<td><a href="https://cran.r-project.org/web/packages/UpSetR/index.html">VennDiagram,UpSetR</a></td>
        </tr>
		 <tr>
            <td>GO enrichment analysis</td>
            <td>GO enrichment analysis for TS related genes</td>
            <td>TS events matrix and group ID</td>
            <td>PDF file and GO result</td>	
			<td><a href="https://bioconductor.org/packages/release/bioc/html/topGO.html">topGO</a></td>
        </tr>
		 <tr>
            <td>Alternative splicing</td>
            <td>Alternative splicing for TS transcript pairs</td>
            <td>TS events matrix</td>
            <td>The summary and deatil information of alternative splicing </td>	
			<td><a href="https://bioconductor.org/packages/release/bioc/html/IsoformSwitchAnalyzeR.html">IsoformSwitchAnalyzeR</a></td>
        </tr>
		 <tr>
            <td>Single gene TS events Identification</td>
            <td>identify TS events for time-series and population transcriptome comparison in single gene level</td>
            <td>sQlite file and gene ID </td>
            <td>html report</td>	
			<td><a href="https://github.com/wyguo/TSIS">TSIS</a>,<a href="http://www.zzlab.net/GAPIT/">GAPIT</a></td>
        </tr>
        <tbody>
    </table>

