# Antarctica_scripts
Scripts used as part of metagenome/16S data analysis

Contained within this repository are scripts used for the following procedures:

## CD-HIT read count clustering
### Read count clustering
*cluster_read_counts.py*
This program reads in a cdhit .clstr file and htseq output file and clusters reads from redundant mappings into a single count corresponding to the representative sequence from each cluster. This script is intended as a precursor to the clustered_read_counts_combine.py script which will take individually clustered files and create a table that edgeR can read. This script is likely to be a bit heavy on memory depending on file size. 

*clustered_read_counts_combine.py*
This program reads in a group of clustered read counts files originating from the cluster_read_counts.py script and converts this into a table that is readable by edgeR. As with the previous script, this could be memory intensive depending on file sizes.

## goseq related
### Signature terms from goseq
*goseq_signature_terms.py*
This script analyses a directory containing tab-delimited text files of all pairwise comparisons of treatments. These files must contain two columns, where column 1 contains functional terms, and column 2 contains the associated FDR/P-values for these. The file's header will be ignored if it is present. The script is also able to have a FDR/P-value cut-off enforced when obtaining signature terms.
