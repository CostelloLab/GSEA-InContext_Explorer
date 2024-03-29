# GSEA-InContext Explorer: the Shiny app for gene set enrichment analysis using GSEA-InContext

Visit the web app [here](http://gsea-incontext_explorer.ngrok.io/) 

GSEA-InContext Explorer allows users to perform two methods of gene set enrichment analysis (GSEA). The first, GSEAPreranked, applies the GSEA algorithm in which statistical significance is estimated from a null distribution of enrichment scores generated for randomly permuted gene sets. The second, GSEA-InContext, incorporates a user-defined set of background experiments to define the null distribution and calculate statistical significance. GSEA-InContext Explorer allows the user to build custom background sets from a compendium of over 5,700 curated experiments, run both GSEAPreranked and GSEA-InContext on their own uploaded experiment, and explore the results using an interactive interface. This tool will allow researchers to visualize gene sets that are commonly enriched across experiments and identify gene sets that are uniquely significant in their experiment, thus complementing current methods for interpreting gene set enrichment results.

Note: Some data files were too large to include in this repo. They can be found in Synapse here: https://www.synapse.org/#!Synapse:syn18898488/files/
