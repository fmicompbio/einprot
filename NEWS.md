# einprot 0.6.4

* Include test results in returned SummarizedExperiment

# einprot 0.6.3

* Combine PTM tests into a single function

# einprot 0.6.2

* Lift out PTM tests into separate functions
* Change interface to defining protein ID column in PTM workflow
* Add tour to iSEE script

# einprot 0.6.1

* Use only spike features with no missing values for normalization
* Add stringIdCol argument
* Lift out filtering and SA plots into separate functions

# einprot 0.6.0

* Change handling of feature identifiers - new arguments idCol, labelCol

# einprot 0.5.22

* Check STRING connection, skip STRING analysis in case of problems

# einprot 0.5.21

* Add idCol argument to PD-TMT workflow

# einprot 0.5.20

* Add samSignificance argument to runTest

# einprot 0.5.19

* Bugfix in PTM workflow to correctly access comparisons
* Add function to map UniProtIDs to gene symbols
* Include STRING networks also for PD-TMT workflow
* Allow t-test for PD-TMT workflow

# einprot 0.5.18

* Allow also Abundances.grouped for PD data

# einprot 0.5.17

* Allow custom names of comparisons

# einprot 0.5.16

* Show only features that are significant in at least one comparison in upset plot
* Expand capabilities of test (a group can belong to more than one merged group, a group can be tested against its complement)

# einprot 0.5.15

* Move normalization to a separate function
* Allow estimating normalization factors from spike features
* Move normalization before imputation

# einprot 0.5.14

* Use eBayes rather than treat if minlFC = 0
* Allow using proDA for differential abundance testing
* Allow the user to specify which PD quantification file to use
* Allow running PD workflow without pdAnalysis file
* Add first version of PD-TMT PTM workflow

# einprot 0.5.13

* Update readFragPipeInfo to use log file instead of config file

# einprot 0.5.12

* Swap UpsetR for ComplexUpset for simplicity
* Include individual volcano plots in table of contents

# einprot 0.5.11

* Update link to download CORUM complexes

# einprot 0.5.10

* Reset paths
* Switch to ComplexUpset for filtering plot
* Use PMID for custom complexes if available

# einprot 0.5.9

* Add decoy filtering for FragPipe output

# einprot 0.5.8

* Don't do upset plots with only one non-empty set

# einprot 0.5.7

* Add FragPipe analysis workflow

# einprot 0.5.6

* Allow subtraction of background as a means of batch correction

# einprot 0.5.5

* Expand description of EMM plots
* Move references to bib file

# einprot 0.5.4

* Add support for reading FragPipe quantifications

# einprot 0.5.3

* Allow specifying the columns to use as primary/secondary feature IDs

# einprot 0.5.2

* Allow running MQ workflow without XML file

# einprot 0.5.1

* Bugfix for features with p-value = 0 in the volcano plots

# einprot 0.5.0

* Switch to a more modular implementation

# einprot 0.4.6

* Allow providing a sample annotation table in TMT workflow

# einprot 0.4.5

* Allow providing a sample annotation table in MQ workflow

# einprot 0.4.4

* Display (up to 6) samples with different symbols in complexes barplots

# einprot 0.4.3

* Bugfix - use correct filter column names in TMT workflow

# einprot 0.4.2

* Allow filtering by (Sequest HT) score also for TMT workflow
* Automatically generate QC pdf for TMT workflow
* Make PCA plots interactive
* Export volcano plots for significant complexes in TMT workflow
* Add possibility to include interactive volcano plots in MQ workflow
* Add direct links to PomBase and WormBase
* Add link table to TMT workflow

# einprot 0.4.1

* Extract more information from the pdAnalysis file
* Add possibility to include interactive volcano plots in TMT workflow
* Add hierarchical clustering of samples to TMT workflow
* Export all feature collection test results to text files
* Subset all feature collections to features in the filtered data set (MQ + TMT)

# einprot 0.4.0

* Add PD-TMT workflow

# einprot 0.3.6

* Export text file with iBAQ values (+mean/sd) for all significant features (across comparisons)
* Add table with direct links to UniProt pages for majority protein IDs

# einprot 0.3.5

* Update STRING to v11.5

# einprot 0.3.4

* Add function to list available complex DBs

# einprot 0.3.3

* Add STRING plots for each comparison
* Shorten column labels in exported heatmap

# einprot 0.3.2

* Export centered heatmap to pdf
* Add upset plot for significant proteins across comparisons
* Add mergeGroups argument to create merged sample groups

# einprot 0.3.1

* Use semicolon instead of comma as separator for PMIDs in complex db
* Include additional complex info in camera output also for t-test

# einprot 0.3.0

* Add more checks of input arguments to runMaxQuantAnalysis
* Expand feature collections to include all variants found in the data set for a given gene
* Internal improvements for increased robustness

# einprot 0.2.2

* Fix display of samples to include/exclude in summary table

# einprot 0.2.1

* Bugfix: fix parsing error for long lines

# einprot 0.2.0

* Collapse duplicates in complex db
* Add path to complex db as argument to runMaxQuantAnalysis

# einprot 0.1.1

* Change path to submission list
* Rotate axis labels in complex bar plots

# einprot 0.1.0

* Initial version
