import gsea_incontext
import sys, logging

def test_function(string):
	return(string)

rnk = "gsea_incontext_test_data/GSE5145_DEG_Expt1_Control_vs_Group1_gene.rnk"
gene_sets = "gsea_incontext_test_data/hallmarks.gmt"
outdir = "out"

# Run Preranked
def run_gsea_preranked(rnk, gene_sets, outdir, permutation_num=100, processes=4):
	preranked_results = gsea_incontext.prerank(
		rnk=rnk,
		gene_sets=gene_sets,
		outdir=outdir, 
		permutation_num=permutation_num,
		processes=processes
	)
	return(True)

run_gsea_preranked(rnk, gene_sets, outdir)