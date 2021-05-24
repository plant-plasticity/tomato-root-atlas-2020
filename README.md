# Tomato-root-atlas-2020
Code and resources for Tomato Root Atlas

**Title:** Novelty, conservation, and repurposing of gene function in root cell type development

**Authors:**  Kaisa Kajala, Mona Gouran, Lidor Shaar-Moshe, G. Alex Mason, Joel Rodriguez-Medina, Dorota Kawa, Germain Pauluzzi, Mauricio Reynoso, Alex Canto-Pastor, Concepcion Manzano, Vincent Lau, Mariana A. S. Artur, Donnelly A. West, Sharon B. Gray, Alexander T. Borowsky, Bryshal P. Moore, Andrew I. Yao, Kevin W. Morimoto, Marko Bajic, Elide Formentin, Niba Nirmal, Alan Rodriguez, Asher Pasha, Roger B. Deal, Daniel Kliebenstein, Torgeir R. Hvidsten, Nicholas J. Provart, Neelima Sinha, Daniel E. Runcie, Julia Bailey-Serres, Siobhan M. Brady

**Journal:** [Cell](https://authors.elsevier.com/sd/article/S0092-8674(21)00504-3)

**Abstract:** Plant species have evolved myriads of solutions, including complex cell type development and regulation, to adapt to dynamic environments. To understand this cellular diversity, we profiled tomato root cell type translatomes. Using xylem differentiation in tomato, examples of functional innovation, repurposing, and conservation of transcription factors are described, relative to the model plant Arabidopsis. Repurposing and innovation of genes are further observed within an exodermis regulatory network and illustrate its function. Comparative translatome analyses of rice, tomato, and Arabidopsis cell populations suggest increased expression conservation of root meristems compared with other homologous populations. In addition, the functions of constitutively expressed genes are more conserved than those of cell type/tissue-enriched genes. These observations suggest that higher order properties of cell type and pan-cell type regulation are evolutionarily conserved between plants and animals.


## Code used by topic (*in progress*)
- [Obtaining RNA-seq read counts: Trimming, Kallisto and STAR](https://github.com/plant-plasticity/tomato-root-atlas-2020/blob/master/Scripts/Trimming_Kallist_STAR_mapping.sh)
- [Relative differential expression: Roku](https://github.com/plant-plasticity/tomato-root-atlas-2020/tree/master/Scripts/Roku)
- [Co-expression network analysis: WGCNA ](https://github.com/plant-plasticity/tomato-root-atlas-2020/blob/master/Scripts/WGCNA-Tomato-ATLAS.rmd)
- [ATAC-aseq analysis](https://github.com/plant-plasticity/tomato-root-atlas-2020/tree/master/Scripts/ATAC)
- [Motif enrichment](https://github.com/plant-plasticity/tomato-root-atlas-2020/tree/master/Scripts/Motif_enrichment)
- [Ontology enrichment](https://github.com/plant-plasticity/tomato-root-atlas-2020/tree/master/Scripts/Ontology_enrichment)
- [Orthology maps](https://github.com/plant-plasticity/tomato-root-atlas-2020/tree/master/Scripts/Orthology_maps)
- [Detecting constitutively expressed genes](https://github.com/plant-plasticity/tomato-root-atlas-2020/tree/master/Scripts/ConstitExpr_homologCT.R)
- [Differential expression: limma-voom](https://github.com/plant-plasticity/tomato-root-atlas-2020/tree/master/Scripts/limmaVoom_DEgenes.R)

## Protocols
These detailed bench protocols describe the methods as carried out in the Tomato Root Atlas paper.
- [TRAP purification of translating ribosomes](https://github.com/plant-plasticity/tomato-root-atlas-2020/blob/master/Protocols/TRAP%20Atlas.pdf)
- [RNA-seq library preparation from TRAP](https://github.com/plant-plasticity/tomato-root-atlas-2020/blob/master/Protocols/Brads%20Rapid%20Ravi%20for%20low%20mRNA%20Atlas.pdf)
- [INTACT isolation of nuclei](https://github.com/plant-plasticity/tomato-root-atlas-2020/blob/master/Protocols/INTACT%20Atlas.pdf)
- [ATAC-seq library preparation from INTACT](https://github.com/plant-plasticity/tomato-root-atlas-2020/blob/master/Protocols/ATAC-seq%20library%20prep%20Atlas.pdf)

## Supplementary data
- [Pairwise comparison of ATAC cut-counts between cell types](https://github.com/plant-plasticity/tomato-root-atlas-2020/blob/master/Figures/Figure_S23_scatter_plot_replicates_repUnion_THSs_ALL_110520_v2_with_legend.pdf)

## Other resources by topic
- [Data on GEO](http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE149217)
- [eFP browser](http://bar.utoronto.ca/eplant_tomato/) Click on the “Tissue and Experiment eFP Viewers.”


## Analysis pipeline summaries

- **Flowchart for identification and annotation of ATAC-seq data.** Data analysis overview for methods used for identification of replicate transposase hypersensitive sites. 

   <img src="https://github.com/plant-plasticity/tomato-root-atlas-2020/blob/master/Figures/Figure_S20_ATACseq_flowchart.jpg" width="700">


- **Unique cell type network construction pipeline.** Data analysis overview for methods used to construct inferred cell type-unique regulatory networks. 

   <img src="https://github.com/plant-plasticity/tomato-root-atlas-2020/blob/master/Figures/Figure_S26_NetworkFlowchart.jpg" width="700">

- **Pipeline of cross-species analysis.** An overview of the pipeline and methods used for exploring the conservation of homologous cell types and tissues among tomato, Arabidopsis, and rice.

  <img src="https://github.com/plant-plasticity/tomato-root-atlas-2020/blob/master/Figures/Cross_species_analysis_overview.jpg" width="700">
  
