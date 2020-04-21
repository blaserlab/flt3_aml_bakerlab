# flt3_aml_bakerlab

This is the analysis repository for 10X data from Daelynn and the Baker lab relating to 4 Flt3+ AML patients treated with gilteritinib.  

See report_20200326.pdf for a summary of the analysis.

High resolution (pdf) plots are available in the directory plots_out/.

Summary and top gene stats are in the directory data_out/.

Feel free to comment by raising an issue, making a pull request, commenting on a commit or sending me an email.

--------------------------------------
## Blog

### April 11 2020

I recalculated the cell counts in each partition after normalizing for the number of cells recovered from each sample.  These can be found here:  "data_out/sample_pa_counts.csv".  This replaces the old file "partition_assignment_counts.csv".  In the new file you want to look at the columns "norm_cell_num" and "norm_pct".

I remade the gene module heatmaps to make the column colors easier to tell apart.  See "plots_out/modules_granular_noTBNKPC.pdf".  I also remade the heatmap including only the AML/Progenitor and Blast 1 partitions.  See "plots_out/modlues_aml_blast1.pdf".  Let me know if you have any questions about these.

I wasn't able to find a 10X dataset for normal bone marrow to use as reference.  I'll keep looking and hopefully have something by the time we get another batch of samples.

Let me know if you have any questions!

Brad

### April 21 2020
#### Question:  
I was wondering if there is a way to just look at the top genes that changed between pre-Gilt treatment and post-Gilt treatment for each patient?  

#### Response:  
I took the main dataset we have been working with. I got rid of T, B, NK, Plasma cells.  Then I stratified everything into 2 timepoint bins:  pre-treatment and post-treatment.  I then calculated the top 50 specific genes for each patient pre- and post-treatment.  They are in the directory "data_out/topgenes_prepost".  Let me know if this isn't what you wanted.

#### Question:  
(Regarding gene modules)...but would we only see over big changes in multiple genes rather than individual changes in genes that may be larger?

#### Response:  
I'm not sure I quite understand - I think there may have been a typo there.  But if your question is "Might individual genes within a module be more or less variable between cells?"  The answer is yes.  To get the overall score for the module, which corresponds to the color on the heatmap, the computer first identifies the gene modules, then for each cell adds up a scaled value for each module, the groups the cells in whatever way we tell it to, and then averages the scaled, summed module score for each group of cells.

The nice thing about this is that if you pick meaningful cell groupings, you can expect there to be differences in module scores between the groups.  Imagine if you divided the dataset randomly into, say, 10 groups of cells.  Then you are averaging together module scores for AML cells, lymphocytes...etc.  And all of the modules end up the same.  However if your groupings have biological meaning, then the average score in these groups of cells will be different.  

The bad thing about modules is that it is an average composite score.  So if you pick cells from two different groups, on average you would expect the module score to predict the expression of module genes in those cells.  However you can't say for certain that any particular gene is up or down in those cells based on the module score.  So it is complementary to the other methods we use to identify top markers.  The bigger the difference in module score between two groups, the more likely any individual gene in that module is to be different between those two groups, if that makes sense.  
