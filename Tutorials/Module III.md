# **Module III: Exploration of TS events**
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
            <td>Filter TS events</td>
            <td>Inital TS matrix</td>
            <td>Filtered TS matrix</td>
			<td>In-house scripts</td>
        </tr>
        <tr>
            <td>Feature-driven exploration</td>
            <td>Characterizing TS events at the genome, gene, transcript, and protein levels</td>
            <td>TS matrix and feature matrix</td>
            <td>HTML report</td>	
			<td>In-house scripts</td>
        </tr>
		 <tr>
            <td>Visualization of expression patterns</td>
            <td>Display expression pattern, exon-intron structure and protein domain of paired transcripts in TS events</td>
            <td>TS matrix (sQlite file), and transcript ID </td>
            <td>PDF file</td>	
			<td>In-house scripts</td>
        </tr>
		 <tr>
            <td>Venn diagram</td>
            <td>Display similatities, differences and enrichment for different group TS events matrix</td>
            <td>TS matrix</td>
            <td>PDF file</td>	
			<td><a href="https://cran.r-project.org/web/packages/UpSetR/index.html">VennDiagram,UpSetR</a></td>
        </tr>
		 <tr>
            <td>GO enrichment analysis</td>
            <td>GO enrichment analysis for TS-related genes</td>
            <td>TS matrix and group ID</td>
            <td>PDF file and GO result</td>	
			<td><a href="https://bioconductor.org/packages/release/bioc/html/topGO.html">topGO</a></td>
        </tr>
		 <tr>
            <td>Alternative splicing</td>
            <td>Alternative splicing analysis for paired transcripts in TS events</td>
            <td>TS matrix</td>
            <td>Information of alternative splicing events</td>	
			<td><a href="https://bioconductor.org/packages/release/bioc/html/IsoformSwitchAnalyzeR.html">IsoformSwitchAnalyzeR</a></td>
        </tr>
		 <tr>
            <td>Identification of TS events for a single gene</td>
            <td>Identify TS events for a user-specified gene</td>
            <td>sQlite file and gene ID</td>
            <td>HTML report</td>	
			<td><a href="https://github.com/wyguo/TSIS">TSIS</a>,<a href="http://www.zzlab.net/GAPIT/">GAPIT</a></td>
        </tr>
        <tbody>
 </table>