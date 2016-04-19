# Tongue vs tongue sample migration plots

Run using QIIME:

```
beta_diversity_through_plots.py -i otu_t_A.biom -m metadata_tongue_vs_tongue_30.txt -o Rare1_659 -t fasttree_all_seed_OTUs.tre
beta_diversity_through_plots.py -i otu_t_B.biom -m metadata_tongue_vs_tongue_30.txt -o Rare2_659 -t fasttree_all_seed_OTUs.tre
```

Run in terminal:

```
Rscript --vanilla plot_differences.R
```

The plot is output into `UniFrac_tvst_movement.pdf`

