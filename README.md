OCP Tools
--
Set of OCP and Focus Panel tools to process NCI-MATCH data.  This set of tools requires several packages from other repos and sources to be installed.
Will update docs with more requirements and maybe even include an installer at some point! 

Current set of scripts and programs available:

   * **ocp_cnv_report.pl**:         Generate a CNV report from a VCF file containing IR CNV data.  Can filter by gene or CN amplitude.

   * **ocp_fusion_report.pl**:      Generate a report of fusions detected in an OCP VCF file.  Can show data for whole panel or just positives.

   * **ocp_control_summary.pl**:    Generate a summary report of the expression control reads in an OCP VCF file.

   * **match_moi_report.pl**:       Run rules to generate a report of NCI-MATCH MOIs and for a VCF file.

   * **variant_review.py**:         Python wrapper script to generate a variant review analysis directory starting with a DNA and an RNA BAM file.

   * **match_control_report.pl**:   Input one or more VCF files from a MATCH control run and output a report.

   * **collate_moi_reports.py**:    Concatenate a group of MOI reports for comparison analysis downstream.  

   * **get_metrics_from_vcf.py**:   Get some quality metrics from VCF or set of VCFs for reporting.

   * **match_delinker.py**:         Script to delink MATCH data for use in other studies.  

   * **ocp_gene_expression_summary.pl**:  Generate some gene expression data from OCPv3 fusion panel assay.  Will allow for better RNA panel quality checks.
   
See the help documentation for each (`<program_name> --help`) for more detailed information about each. 
