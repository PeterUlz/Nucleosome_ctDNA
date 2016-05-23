# Inferring expressed genes by whole-genome sequencing of plasma DNA

The analysis of cell-free DNA (cfDNA) in plasma represents a rapidly 
advancing field in medicine. cfDNA consists predominantly of nucleosome-protected
 DNA shed into the bloodstream by cells undergoing apoptosis. We performed whole-genome
 sequencing (WGS) of plasma DNA and identified two discrete regions at transcription 
start sites (TSS) where the nucleosome occupancy results in different read-depth coverage
 patterns in expressed and silent genes. By employing machine learning for gene classification,
 we found that the plasma DNA read depth patterns from healthy donors reflected the expression 
signature of hematopoietic cells. In cancer patients with metastatic disease, we were able to
 classify expressed cancer driver genes in regions with somatic copy number gains with high accuracy.
 We could even determine the expressed isoform of genes with several TSSs as confirmed by RNA-Seq 
of the matching primary tumor. 
Our analyses provide functional information about the cells releasing their DNA into the circulation.  

A preprint can be found on [bioRxiv](http://biorxiv.org/content/early/2016/04/20/049478/ "bioRxiv Preprint")  

This is a collection of scripts which were used to analyse whole-genome sequencing data of non-cancer controls
as well as two breast cancer samples.

* Alignment and read trimming
* Coverage at chromosome 12 nucleosome array
* Average coverage around Transcription start sites
* Expression prediction based on two features derived from the average coverage analysis
    * Coverage from -1000bp to +1000bp from TSS (2K-TSS Coverage)
    * Coverage from -150bp to +50bp from TSS (NDR Coverage)


