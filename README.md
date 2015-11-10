OCP Tools
--
Set of OCP and Focus Panel tools to process NCI-MATCH data.  This set of tools requires several packages to be installed, the 
most important being Thermo Fisher's VCF converter script. Current set of scripts and programs available. 

   * **oncomine-vcf-converter**:    Thermo Fisher script to parse out a VCF file into a TSV file that can be more easily consumed downstream.
   * **ocp_cnv_report.pl**:         Generate a CNV report from a VCF file containing IR CNV data.  Can filter by gene or CN amplitude.
   * **ocp_fusion_report.pl**:      Generate a report of fusions detected in an OCP VCF file.  Can show data for whole panel or just positives.
   * **ocp_control_summary.pl**:    Generate a summary report of the expression control reads in an OCP VCF file.
   * **match_moi_report.pl**:       Run variant rules to generate a report of NCI-MATCH MOIs and aMOIs for a VCF file.  
   * **plotCNV.R**:                 Thermo Fisher R script to generate a plot of CNV results for an OCP panel.
   * **variant_review.py**:         Python wrapper script to generate a variant review analysis directory starting with a DNA and an RNA BAM file.
   
See the help documentation for each (`<program_name> --help`) for more detailed information about each. 
