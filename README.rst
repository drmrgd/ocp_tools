#########
OCP Tools
#########

************
Introduction
************

Set of Oncomine Cancer Panel OCP (renamed as Oncomine Comprehensive Assay (OCA))
and Focus Panel tools to process NCI-MATCH data.  This set of tools requires 
several packages from other repos and sources to be installed.

This set of scripts requires the use of the `ir_utils package <https://github.com/drmrgd/ir_utils.git>`_ as well as samtools and other common libraries and programs.
More script level requirements listed below with each script.

*********************************************
Current set of scripts and programs available
*********************************************

   * **collate_moi_reports.py**:
       - Concatenate a group of MOI reports generated with ``match_moi_report.pl``
         for comparison analysis downstream. A bit primitive, but can be helpful
         for quickie large analyses.

   * **get_metrics_from_vcf.py**:
       - Get some quality metrics from VCF or set of VCFs for reporting.  Can 
         report on MAPD, RNA reads, and expression control data.

   * **match_delinker.py**:
       - Script to delink MATCH data for use in other studies.

   * **match_moi_report.pl**:
       - Run rules to generate a report of NCI-MATCH MOIs for a NCI-MATCH VCF file.

   * **match_positive_control_report.pl**:
       - Input one or more VCF files from a MATCH control run and output a report.

   * **match_rna_qc.pl**:
       - New to panel is multi-pooled RNA assays.  This script will read the pool
         level controls to determine how the panel performed as a whole. Requires 
         proprietary panel JSON file that is not for distribution. You can 
         generate this yourself by creating a JSON file starting from the 
         fusions BED file that accompanies each assay panel, and includes the
         following stucture: ::

           {
               "ExprControl" : {
                   "assay_id" : "pool_number"
               },
               "GeneExpression" : {
                   "assay_id" : "pool_number"
               },
               "Fusion" : {
                   "assay_id" : "pool_number"
               }
           }

   * **ocp_cnv_report.pl**:
       - Generate a CNV report from a VCF file containing IR CNV data.  Can 
         filter by gene or CN amplitude. One component of ``match_moi_report.pl``.

   * **ocp_control_summary.pl**:
       - Generate a summary report of the expression control reads in an OCP 
         VCF file.  Deprecated and will be replaced fully at some point by the
         ``match_rna_qc.pl`` script.

   * **ocp_fusion_report.pl**:
       - Generate a report of fusions detected in an OCP VCF file.  Can show 
         data for whole panel or just positives.

   * **variant_review.py**:
       - Python wrapper script to generate a variant review analysis directory 
         starting with a DNA and an RNA BAM file.  This wrapper requires the 
         ``ir_utils`` pacakge in order to run.

   * **dl_reporter.py**:
       - Python3 script that relies on the ``matchbox_api_utils`` package to 
         add aMOI annotations to a MOI report.  This is useful for the so-called
         Designated Labs effort, where the upfront data is not going to be run
         through MATCHBox for matching.  Right now, relies on 
         ``match_moi_report.pl`` (which is why it's a part of this repo), but 
         ultimately should be capable of reading any data.

See the help documentation for each (``<program_name> --help``) for more detailed
information about each. 
