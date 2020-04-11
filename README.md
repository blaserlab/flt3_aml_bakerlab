# flt3_aml_bakerlab

This is the analysis repository for 10X data from Daelynn and the Baker lab relating to 4 Flt3+ AML patients treated with gilteritinib.  

See report_20200326.pdf for a summary of the analysis.

High resolution (pdf) plots are available in the directory plots_out/.

Summary and top gene stats are in the directory data_out/.

Feel free to comment by raising an issue, making a pull request, commenting on a commit or sending me an email.

--------------------------------------

## Update April 11 2020

I recalculated the cell counts in each partition after normalizing for the number of cells recovered from each sample.  These can be found here:  "data_out/sample_pa_counts.csv".  This replaces the old file "partition_assignment_counts.csv".  In the new file you want to look at the columns "norm_cell_num" and "norm_pct".

I remade the gene module heatmaps to make the column colors easier to tell apart.  See "plots_out/modules_granular_noTBNKPC.pdf".  I also remade the heatmap including only the AML/Progenitor and Blast 1 partitions.  See "plots_out/modlues_aml_blast1.pdf".  Let me know if you have any questions about these.

I wasn't able to find a 10X dataset for normal bone marrow to use as reference.  I'll keep looking and hopefully have something by the time we get another batch of samples.

Let me know if you have any questions!

Brad
