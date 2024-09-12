# This is the development version of the R package _phenotools_

Please report any bugs or issues, as this package is currently under development.

## Concept

Here is a 

![Figure 2-01](fig2-01.jpg)

Character similarity networks and character annotation. Anatomy of a character statement (A). Arrows indicate two approaches for comparing characters: fuzzy text matching of terms (lines without arrowheads) and text searching within characters (lines with arrowheads). Network diagrams in gray boxes show connections among characters using different approaches, with thicker lines indicating more similar characters. (B) Automated semantic character annotation with generate_ontology. Character statements can be structured along a continuum, with semi-structured characters (Clarke, Zhou, & Zhang, 2006) presenting greater opportunity for automated phenomic dataset assembly.

## Installing _phenotools_

You can install the package using the devtools R package:
```
install.packages("devtools")
library(devtools)
install_github("celiason/phenotools")
library(phenotools)
```
