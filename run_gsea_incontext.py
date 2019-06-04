import gsea_incontext
import sys, logging

def test_function(string):
	return(string)

#rnk = 'gsea_incontext_test_data/GSE5145_DEG_Expt1_Control_vs_Group1_gene.rnk'
#gene_sets = 'data/gene_sets/hallmarks.gmt'
#background_csv = 'out/incontext_bg_input.csv'
#outdir = 'out'

# Run GSEA-InContext
def run_gsea_incontext(rnk, gene_sets, background_csv, outdir, permutation_num=100, no_plot=True, processes=4):
	gseaincontext_results = gsea_incontext.incontext(
		rnk=rnk,
		gene_sets=gene_sets, 
		background_rnks=background_csv, 
		outdir=outdir, 
		permutation_num=permutation_num,
		no_plot=no_plot,
		processes=processes
	)
	return(True)

#run_gsea_incontext(rnk, gene_sets, background_csv, outdir)