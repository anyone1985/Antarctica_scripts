# Antarctica_scripts
Scripts used as part of metagenome/16S data analysis

Contained within this repository are scripts used for the following procedures:

## CD-HIT related
### Read count clustering
*cluster_read_counts.py*

This program reads in a cdhit .clstr file and htseq output file and clusters reads from redundant mappings into a single count corresponding to the representative sequence from each cluster. This script is intended as a precursor to the *clustered_read_counts_combine.py* script which will take individually clustered files and create a table that edgeR can read. This script is likely to be a bit heavy on memory depending on file size. 

*clustered_read_counts_combine.py*

This program reads in a group of clustered read counts files originating from the *cluster_read_counts.py* script and converts this into a table that is readable by edgeR. As with the previous script, this could be memory intensive depending on file sizes.

## GOseq related
### Signature terms from GOseq
*goseq_signature_terms.py*

This script analyses a directory containing tab-delimited text files of all pairwise comparisons of treatments. These files must contain two columns, where column 1 contains functional terms, and column 2 contains the associated FDR/P-values for these. The file's header will be ignored if it is present. The script is also able to have a FDR/P-value cut-off enforced when obtaining signature terms.

## QIIME related
### NMDS table creation
*qiime_table_to_nmds.py*

This program will read in the table files QIIME creates from the taxa_summary_plots script and reformat this data into a table that can be read into R and VEGAN.

### Krona file formatting
*qiime_to_ktImportText.py*

This program will read in the table files QIIME creates from the taxa_summary_plots script and reformat this data into a table that can be read by ktImportText to produce Krona plots.

## k-SLAM related
### File preparation for downstream analysis
*kslam_xml_taxIdMap.py*

This program will read in the xml files k-SLAM produces and create a tab-delimited file that can be used for later Krona or NMDS plot generation. Note that the main benefit of producing a Krona plot through this avenue (rather than by using ktImportTaxonomy) is that file size can be significantly reduced. The plot is not as well-formatted as the one ktImportTaxonomy produces by default, however.

*kslam_taxIdMap_to_ktImportText.py*

Reads in the output(s) of *kslam_xml_taxIdMap.py* to produce tab-delimited file(s) that can be called by ktImportText to produce Krona plots.

*kslam_taxIdMap_to_nmds_group.py*

Reads in the output(s) of *kslam_xml_taxIdMap.py* to produce a tab-delimited file that can be read into R and used by VEGAN to produce NMDS plots. The taxonomic rank to summarise abundances to can be specified if a NMDS at the deepest rank assigned by k-SLAM is not wanted.
