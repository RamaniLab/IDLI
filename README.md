# Iteratively Defined Lengths of Inaccessibility (IDLI):

This README is intended to share code for analyses performed in the following manuscript:

Yang MG, Wang S, Richter HJ, McNally CP, Harris N, Maresca M, de Wit E, Willenbring H, Maher J, Goodarzi H#, and Ramani V#. "Pervasive, patterned, and programmed nucleosomal distortion on mammalian chromatin fibers". *Submitted* (2024)

*All scripts are provided as is, without any warranty and for use at your own risk. This is __not__ the release of a software package. We are only providing this information and code, in addition to a description of methods, for making it easier to reproduce our analyses. We are __not__ providing any support for these scripts.*

Summary:
--------------------

[ABSTRACT HERE]

Code overview:
--------------------
*Preprocessing PacBio data for single-molecule accessibility*

To preprocess PacBio data, we utilized software from Pacific Biosciences and a custom script ([hmm_output_t_values.py](https://github.com/RamaniLab/Subnucleosome/blob/main/hmm_output_t_values.py)) to successively run our NN-HMM across a range of user-specified transition probabilities. These resulted in the following output per sample: (1) an alignment of CCS reads to the mouse reference genome (mm10); (2) an accessibility prediction per CCS molecule (Viterbi path of HMM component of model); and (3) an atlas of footprint locations (*i.e.* start / end positions) and sizes on a per-molecule basis. Files for (2) and (3) were generated for *t*  = 1, 31, 51, 71, and 101 / 1,000 for each sequencing library.

*Quantification of extent of methyltransferase footprinting per sequencing library*

We utilized the ratio of mononucleosome-sized to dinucleosome-sized footprints (MDR) as a proxy for measuring the extent of EcoGII footprinting on a per-sample basis. ‘Poorly methylated’ samples exhibit a higher fraction of footprints longer than the expected length of DNA protected by a mononucleosome, which can be ascribed to the failure of EcoGII to efficiently methylate in linker DNA. We employed a custom script ([compute_mdr_per_sample.py](https://github.com/RamaniLab/Subnucleosome/blob/main/compute_mdr_per_sample.py)), which uses the *find_peaks* function (from *scipy*) and computes the ratio of maximum peak heights in footprint length histograms for mononucleosome- and dinucleosome-sized footprints per sequencing library.

*Visualization of single-molecule data across different t values at single genomic loci*

All datasets from E14 mESCs with high MDR (>10) were processed with a custom script ([generate_bam_for_igv.py](https://github.com/RamaniLab/Subnucleosome/blob/main/generate_bam_for_igv.py)) that encodes single-molecule accessibility patterns (for any given *t* value) within a .bam file with alignment information. SAMOSA data at the SOX2 locus (chr3: 34,756,929 – 34,759,161; including the SOX2 promoter, gene body, and downstream SOX2 Control Region (SCR)) were visualized with a custom version of IGV, as previously described (https://github.com/RamaniLab/SMRT-Tag/tree/main/igv-vis). Bulk ATAC-seq data from mESCs were obtained as processed .bw files (GSE98390) and visualized at identical coordinates in the UCSC Genome Browser.

*Clustering single-molecule accessibility patterns with respect to footprint midpoints*

To examine patterns of nucleosomal distortion, we clustered accessibility patterns across a range of t values with respect to footprint midpoints (i.e. predicted dyad positions). Footprints <200 nt in length (from *t* = 1 / 1,000) within +/- 500 bp of midpoints of epigenomic domains, repeat sequences, or TF-binding motifs were selected. Per-footprint accessibility data (for *t* = 31, 51, and 71 / 1,000) from a 140-bp window centered on each footprint midpoint (*i.e.* 420 bp in total per footprint) were used as input for Leiden clustering (*res* = 0.5, 0.6, 0,6, and 0.6 for domains / repeats, SOX2, CTCF, and FOXA2 motifs, respectively). Clusters were filtered such that the least abundant clusters that collectively summed up to <10% of each dataset were removed. Representative code for performing this analysis is provided ([process_footprints_nucleosomal_distortion.py](https://github.com/RamaniLab/Subnucleosome/blob/main/process_footprints_nucleosomal_distortion.py)).

*Clustering single-molecule accessibility patterns with respect to sequence features of interest*

To define translational settings, we clustered accessibility patterns across a range of t values with respect to midpoints of a given sequence feature of interest (e.g. repeat elements, TF-binding motifs). Footprints <200 nt in length (from *t* = 1 / 1,000) wherein the absolute distance between the footprint midpoint and sequence feature was than half the total length of the footprint in question were selected (e.g. within +/- 100 bp for a 200 nt footprint). Per-footprint accessibility data (for *t* = 31, 51, and 71 / 1,000) from a 140-bp window centered on each sequence feature (*i.e.* 420 bp in total per footprint) were used as input for Leiden clustering (*res* = 0.3, 0.8, 0.8, and 0.8 for repeats, SOX2, CTCF, and FOXA2 motifs, respectively). Clusters were filtered such that the least abundant clusters that collectively summed up to <10% of each dataset were removed. Representative code for performing this analysis is provided ([process_footprints_translational_positions.py](https://github.com/RamaniLab/Subnucleosome/blob/main/process_footprints_translational_positions.py)).

*Footprint size distributions, mean accessibility patterns, and ‘horizon plots’ for subfootprint organization per cluster*

To assign putative biological function to clusters of nucleosomal distortion patterns or translational positions, we visualized the following features. First, we plotted the distribution of footprint sizes on a per-cluster basis for data from the lowest *t* value examined (*t* = 1 / 1,000). Second, we plotted mean accessibility as a function of distance to respective footprint / sequence feature midpoints. These data were specifically selected from the highest *t* value used for Leiden clustering in the previous step (*t* = 71 / 1,000). Third, we generated ‘horizon plots’ (analogous to commonly used ‘V plots’ in MNase-seq data) by computing z-scores for enrichment of footprints of particular size (20 – 200 nt) at particular positions (+/- 100 bp from footprint / sequence feature midpoints) for each individual cluster. z-scores are summed per position across a range of t values (*t* = 31, 51, and 71 / 1,000). Color scale limits for the resultant heatmap are bounded by the 0th (lower) and 95th (upper) percentile z-scores for positions within +/- 70 bp of the footprint / sequence feature midpoint. Instances of clusters with entirely unmethylated nucleosomes were manually filtered out before downstream quantification. Representative code for performing this analysis is provided ([plot_footlen_meanacc_horizondata.py](https://github.com/RamaniLab/Subnucleosome/blob/main/plot_footlen_meanacc_horizondata.py)).

*Structural visualization of nucleosomal accessibility patterns*

For structural visualization of mean accessibility patterns for defined clusters, we utilized a custom script ([accessibility_to_pdb_structure.py](https://github.com/RamaniLab/Subnucleosome/blob/main/accessibility_to_pdb_structure.py)) to convert our mean single-molecule accessibility data per cluster of interest from the highest *t* value used for Leiden clustering (*t* = 71 / 1,000) into a format compatible with visualization via ChimeraX (v1.7.1; *attribute* = percentAcc, *match mode* = 1-to-1, and *recipient* = residues). Accessibility values from residues -76 to +74 were overlaid on PDB 7KBD structure.

*Identification of single-molecule co-occupancy patterns for different footprint groups*

Individual fibers with at least three Leiden-classified footprints of nucleosomal distortion patterns within +/- 500 bp of repeat elements / TF motifs of interest were selected. Fibers with footprints >200 nt and/or with footprints in the least abundant 10% of clusters were excluded from further analysis. Consecutive ‘triplet’ footprints were tabulated by computing midpoints and lengths for: (1) the most proximal footprint upstream of the repeat element / TF-binding motif, (2) the central-most footprint nearest the repeat element / TF-binding motif, and (3) the most proximal footprint downstream of the repeat element / TF-binding motif. Per-molecule accessibility data was computed for a 2-kb window centered on the repeat element / TF-binding motif of interest for visualization alongside footprint positions, accounting for read and sequence feature strand appropriately. Representative code for performing this analysis is provided ([process_triplet_footprints.py](https://github.com/RamaniLab/Subnucleosome/blob/main/process_triplet_footprints.py)).

Data availability:
--------------------
Raw and processed data will be made available at GEO accession GSEXXXXXX.

Contact information:
--------------------
Please contact Vijay Ramani (vijay.ramani[at]gladstone.ucsf.edu) or Marty Yang (marty.yang[at]gladstone.ucsf.edu) for questions.
