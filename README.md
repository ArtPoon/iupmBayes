# iupmBayes

[iupmBayes](http://github.com/ArtPoon/iupmBayes) is a set of R scripts for estimating the number of infectious units per million (IUPM) from HIV sequence data obtained from a quantitative viral outgrowth assay (QVOA).  



## Data

These are two tab-delimited tabular data sets that accompany our manuscript.  They comprise the following variables:

* Participant number
* Well number
* Number of cells plated in the well
* p24 positive (binary outcome)
* Region (HIV gene)
* Total variants for participant by region
* Variants per region per well
* Variant index
* Presence of variant (binary outcome)

These data sets can be read into R using the provided `parse.data` function.
