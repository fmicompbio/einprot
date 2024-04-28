# einprot 0.9.4

* Add details about DIA-NN command line to report
* Add column with links to ComplexPortal query to link table
* Include pathways from MSigDB among supported feature collections
* Update DIA-NN example data
* Expand input data paths in run*Analysis() functions
* Add support for importing Spectronaut PG pivot files

# einprot 0.9.3

* Significantly expand the built-in complex database by adding complexes from the Complex Portal and hu.MAP 2.0
* Hide box plot legends if there are too many groups
* Add MinProbGlobal imputation method
* Add makeInteractiveVolcano argument to plotVolcano

# einprot 0.9.2

* Fix sample column name in DIA-NN import
* Adapt DIA-NN log import to handle vector output
* Add argument addHeatmaps to run functions, to allow suppression of heatmaps for large data sets
* Remove samples with no detected features
* Add support for Spectronaut protein group output
* Add support for displayName column in sample annotation table
* Add extraFeatureCols argument for specification of user-defined feature annotation columns
* Extract information about Top N ions for FragPipe

# einprot 0.9.1

* Replace ggseqlogo dependency with motifStack
* Adapt missing value plots to work with both proportions and percentages

# einprot 0.9.0

* Consolidate tool-specific templates into a single template

# einprot 0.8.0

* Add support for DIA-NN protein group data

# einprot 0.7.7

* Change default value of singleFit to TRUE everywhere, for consistency

# einprot 0.7.6

* Bug fix in filtering plot when only one criterion is available

# einprot 0.7.5

* Move QC plot function

# einprot 0.7.4

* Add argument to allow limiting the labeling to only features with positive/negative logFCs
* Add function to generate heatmaps
* Add function to create summary of abundance values for significant features
* Add argument to define the assays(s) to use for exported values and barplots
* Harmonize treatment of merged groups in runTest and plotVolcano
* Sort exported test results by p-value instead of logFC
* Expand the set of allowed species

# einprot 0.7.3

* Combine download buttons and pagelength selection for tables
* Add check and stop if no features are left after filtering

# einprot 0.7.2

* Add option to export link table to csv
* Add summary table with sequence windows to PTM workflow
* Add shiny app to generate sequence logos for PTMs
* Add link to workflow overview in table of contents
* Include link to top feature sets for each comparison in table of contents

# einprot 0.7.1

* Add scale and subset_row arguments to doPCA
* Move definition of feature IDs after the filtering
* Extract FragPipe decoy pattern from the config file
* Remove features with too few non-imputed values before doing PCA
* Retain a few more rowData columns for MaxQuant
* Fix extraction of fixed modifications for MQ v2.2
* Add max peptide mass and min peptide length to MQ report
* Add Top3 to the valid iColPatterns
* Include Top3 abundances in MQ result tables if available

# einprot 0.7.0

* Change approach for injecting values into the Rmd file, to avoid the need to duplicate escapes
* Allow the user to set the column to colour by in the PCA
* Include 'extra' columns in output from runPTMTest
* Allow '^Abundance.' iColPattern
* Suggest closest iColPattern if the provided one is invalid

# einprot 0.6.10

* Adapt readFragPipeInfo to handle different FragPipe versions
* Add filtering function for FragPipe
* Add option to use sample weights for limma
* Move modifications filters from the PDTMTptm to the PDTMT PeptideGroups workflow
* Add interactiveGroupColumn argument

# einprot 0.6.9

* Safeguard against non-functional connection to STRING server
* Allow disabling score/number of peptides/PTMs filters
* Don't fail if interactiveDisplayColumns are not present in data
* Use Perl-compatible regular expressions when matching link table columns
* Include Sequence column in final SCE for PD-TMT data

# einprot 0.6.8

* Add idCol and labelCol arguments to PTM workflow
* Change default behaviour of fixFeatureIds when column name is given to not make output unique
* Add stringVersion and stringDir arguments, allowing the use of local STRING files
* Make filter functions more robust to missing columns
* Add possibility to write excluded features to a file
* Represent link table columns as integers/factors when appropriate
* Allow displaying any column in rowData(sce) in the interactive volcano plot tooltip
* In the case of long labels, attempt to auto-adapt text size in PCA coefficient and logFC plots
* Add overview and crosslinks in the beginning of the reports
* Include bar plot of significant features in pdf output
* Allow iColPattern without Sample for PDTMT
* Make plot axis labels less ambiguous

# einprot 0.6.7

* Add signifDigits argument to makeDbLinkTable, and round to 4 significant digits in the templates
* Include the maximum number of missed cleavages in PDTMT tables
* Let maxNbrComplexesToPlot determine also the maximum number of top feature sets displayed in the reports
* Bugfix for sample plot ordering in complex bar plots
* Add bar plot for significant features

# einprot 0.6.6

* Bring FragPipe workflows up to date
* Add individual PTM volcano plots to table of content
* Add modificationsCol and keepModifications arguments to PTM workflow
* Increase control level in deparsing to allow e.g. multi-line functions
* Fill feature sets plot by direction
* Suppress legend in PCA plot if there are too many groups
* Add option to only retain master proteins in PDTMT Proteins workflow
* Allow inclusion of extra columns in the link table

# einprot 0.6.5

* Change correlation heatmap appearance
* Provide the possibility to run the workflows without statistical tests
* Allow specifying the iColPattern without escaped periods
* Add interactiveDisplayColumns arguments to volcano plots

# einprot 0.6.4

* Include test results and comparison list in returned SummarizedExperiment
* Expand the PCA plots with scree plots and coefficient plots
* Add option to label any/only significant features in volcano plot
* Add einprotLabel column to link table

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
