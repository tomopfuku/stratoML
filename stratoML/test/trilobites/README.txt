README

Data and Analytical Files for Bapst and Hopkins (pterocephaliids with cal3) manuscript
Prepared 03-25-2016

The supplementary data package for this study consists of several sets of files. The first set of files is the input data files used in the presented analyses.

The pruned most-parsimonious tree of the pterocephaliid trilobites from Hopkins (2011) is provided as a NEXUS file:

	PteroGW-pruned.tre
	
The principal component scores from a PCA of the geometric morphometric landmark data for those pterocephaliid taxa is provided as a tab-delimitted .txt table:

	PCscores-pterocephaliid.txt

The stack of 100 CONOP solutions, consisting of first and last appearance dates for all taxa or for just the pterocephaliids is provided as two tab-delimitted .txt tables:
	
	FAD-LAD-rescaled-all.txt
	FAD-LAD-rescaled-ptero.txt

The minimum number of horizons each taxon was found at (see main text) is also provided as a tab-delimitted .txt tables:
	
	all-minsampling.txt
	ptero-minsampling.txt

All post-inference analysis and visualization were done in R, via an Rmarkdown script in RStudio. The script (a .Rmd), the resulting markdown PDF with output and figures, and the saved workspace file are included in our supplemental data materials:

	trilobite_cal3_03-25-16.Rmd
	trilobite_cal3_03-25-16.pdf
	trilobite_workspace_03-25-2016.Rdata

Finally, all published figures were created entirely from within R, using the saved workspace from the Rmarkdown script and the following R script:
	
	figures_05-30-16.R

All analyses were performed with R v3.2.3 using R packages ape v3.4, geiger v2.06, paleotree v2.6, phangorn v2.02 and mvMorph v1.0.6.