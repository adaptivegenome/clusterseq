##Purpose
ClusterSeq is a tool used to analyze the distribution of barcodes across multiple tagged samples.

##Description

ClusterSeq consists of two separate tools- an executable for filtering and clustering of barcode sequences, and a python script for comparing sequences 

###Cluster command (cluster.cpp)

The cluster command expects reads in a FASTQ format, and expects reads with the following format:

[tag][start_marker][data][end_marker]

Where tag is a multiplexing tag identified in a required taglist file (see format below), and start_marker and end_marker are constant sequences before and after the barcode data. These are specified as command line options.

When run, the cluster command separates each read by tag, then performs the following:

* Check for begin and end 'markers' identifying a valid barcode. Sequences with invalid tags are discarded.
* Cluster barcode sequences with a limited number of differences, in order to tolerate sequencing errors. The number of acceptable differences is specified on the command line.
* Write a file containing a list of all reads and quality strings for this tag. Each line of the file contains the data sequence, a tab, and the corresponding Phred quality.
* Write a file containing a list of barcodes identified, along with the number of similar barcodes found with this tag. This file is a CSV consisting of barcode data followed by the number of times this barcode occurs in the input data (with the current tag).

The taglist file should consist of one or more lines in the following format:

<pre>
[FASTQ_name1] [whitespace] [TAG1]
              [whitespace] [TAG2]
[FASTQ_name2] [whitespace] [TAG1]
              [whitespace] [TAG2]
</pre>

###Cluster comparison tool (merge_clusters.py, merge_clusters2.py)
The cluster merging tool combines the CSV output from multiple tags of the cluster tool so that the frequency of occurrence of barcodes can be compared across multiplexing tags. The output is a CSV file listing the frequency of occurrence of a barcode across each of the input files; the first column lists the barcode, and each subsequent column corresponds to an input file.

When combining clusters from multiple files, you may choose to combine similar clusters. This can be done by withthe cluster_edit_distance_threshold parameter- see below for details.

The size of the generated clusters may also be examined by plotting the data in the merged_clusters_histogram.csv file (see below for details).

The tool generates three files:

* merged_clusters.csv contains all barcodes from all files.
* merged_clusters_filtered.csv contains a subset of the data from merged_clusters.csv to only show barcodes that occur a certain number of times, set by the threshold 'min_count_for_filter' at the top of the script file.
* merged_clusters_histogram.csv contains information needed to generate a histogram of cluster sizes. The first column of this file is the range of sizes included in each bin, and the second column contains the number of clusters containing that number of sequences. This file is designed to be imported into another tool to generate a graph (Excel, Matlab, R, etc).

Two versions of this tool are available, with only a minor difference- 

* 'merge_clusters.py' generates a 'merged_clusters_filtered.csv' file which contains lines in which the barcode appeared at least 'min_count_for_filter' times in **each** input file.
* 'merge_clusters2.py' generates the same file, but with barcodes that appeared at least 'min_count_for_filter' times in **any** input file.

The following parameters can be set at the top of the .py files:

* min_count_for_filter controls filtering of reads. Only reads that occur at least this many times will be included in the output file. Set to zero for no filtering.
* cluster_edit_distance_threshold controls which clusters will be combined. This parameter is the acceptable number of differences between two cluster names; 'N's in either cluster name mean that that position will be ignored for comparison purposes. Set to -1 to not combine clusters.
* histogram_bin_growth_factor controls the size of the bins used when generating the histogram file. Bins thresholds are histogram_bin_growth_factor ^ N, where N is an integer 0 and up. The number of bins depends on the cluster with the largest number of sequences in it, as enough bins will be generated to contain all data.

##Compilation

With OpenMP:
<pre>
g++ cluster.cpp -o cluster -O3 -fopenmp
</pre>

Without OpenMP:
<pre>
g++ cluster.cpp -o cluster -O3
</pre>

##Running
Usage is:
<pre>
cluster min_quality max_n_allowed num_diff_allowed [sample_fraction] FASTQ_name_no_ext [tag_file_name] [start_marker] [end_marker]
</pre>
For example:
<pre>
./cluster 30 2 2 SNP456 taglist.txt AAAAAAAA TTTTTT
</pre>
The name of the FASTQ input file should have no extension- this is also used to load the proper tags of from the taglist file.

Details on the meaning of each of these commands can be found by running the cluster command with no arguments.

The cluster comparison tool simply needs a list files to merge:
<pre>
./merge_clusters.py SNP456.CTAG_clusters.csv SNP456.ACGT_clusters.csv
</pre>

##Author contact
Lee Baker, VBI: leecb@vt.edu

David Mittelman, VBI: david.mittelman@vt.edu
