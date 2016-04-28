# Scripts for generating paper figures

All the figures from the paper Expanding the UniFrac toolbox, submitted to Great Lakes Bioinformatics and the Canadian Computational Biology Conference 2016 (https://www.iscb.org/glbioccbc2016)

The paper can be found here: https://github.com/ruthgrace/unifrac_paper

In this repository, the scripts often use the same data sets, so the data sets are stored separately in the [data folder](data).

## Sample migration in different rarefactions, plotted on principal coordinates, measured with unweighted UniFrac.

[Tongue dorsum sample migration code](tongue_sample_migration)

[Tongue dorsum vs. buccal mucosa sample migration code](tongue_cheek_migration)

Output:

![Sample migration](images_for_README/sample_migration.png)

## Maxiumum relative deviation of rarefactions versus median deviation for traditional and non-traditional microbiome dissimilarity metrics.

[Ellipses plot code](tongue_ellipses)

Output:
![Ellipses plot](images_for_README/ellipses_plot.png)

## Analysis of tongue and buccal mucosa data using different UniFrac weightings.

[PCoA plot code](tongue_dorsum_vs_buccal_mucosa)

Output:

![Tongue dorsum vs buccal mucosa](images_for_README/tongue_cheek_pcoa.png)

## Analysis of breast milk data using different UniFrac weightings.

[PCoA plot code](breastmilk)

Output:

![Breast milk data set with infected sample](images_for_README/breastmilk_pcoa.png)

## Analysis of simulated monocultures using different UniFrac weightings.

[PCoA plot code](monoculture)

Output:

![Monoculture PCoA plots](images_for_README/monoculture_pcoa.png)

## Principal Coordinate Analysis derived from GUniFrac distance matrices.

These figures are in the supporting information section.

[GUniFrac sample migration plot code](supporting_information/gunifrac_migration)

Output:

![GUniFrac PCoA plots](images_for_README/gunifrac_1.png)

![GUniFrac PCoA plots](images_for_README/gunifrac_2.png)

![GUniFrac PCoA plots](images_for_README/gunifrac_3.png)

## Principal coordinate analysis derived from tongue dorsum samples using unweighted UniFrac distance matrices with no tree pruning.

This figure is in the supporting information section.

[Tongue dorsum sample migration without tree pruning plot code](supporting_information/tongue_tongue_no_tree_pruning)

Output:

![Tongue dorsum sample migration without tree pruning](images_for_README/tongue_pruning.png)

## Principal coordinate analysis derived from tongue dorsum and buccal mucosa samples using unweighted UniFrac distance matrices with no tree pruning.

This figure is in the supporting information section.

[Tongue dorsum vs. buccal mucosa sample migration without tree pruning plot code](supporting_information/tongue_cheek_no_tree_pruning)

Output:

![Tongue dorsum vs. buccal mucosa sample migration without tree pruning.](images_for_README/tongue_cheek_pruning.png)
