# Bioinformatics

Bioinformatics is an interdisciplinary field that integrates biology and computer science to effectively analyze and interpret extensive and complex datasets.

# Data types: 
Genomic Data (including WGS and WES):

Whole Genome Sequencing (WGS) Data: Comprehensive sequencing of the entire DNA of an organism.

Whole Exome Sequencing (WES) Data: Focused sequencing of all protein-coding regions in the genome.

Transcriptomic Data: Analysis of RNA transcripts to understand gene expression patterns.

Proteomic Data: Study of the full set of proteins produced in a biological system.

Metabolomic Data: Analysis of metabolites to understand metabolic processes and organism responses.

Epigenetic Data: Study of heritable changes in gene function that do not involve changes in the DNA sequence.

Phylogenetic Data: Analysis of evolutionary relationships among organisms or genes.

Structural Biology Data: Information about the 3D structure of biological molecules.

Microbiome Data: Study of microbial communities in specific environments.

Population Genetics Data: Genetic variation across populations, useful for understanding evolutionary processes.

Clinical Data: Patient-specific data combining genomic, proteomic, and other omics data.

Chemoinformatics Data: Study of chemical substances and their interactions in biological systems.

Pathway Data: Information on biological pathways, key for understanding cellular processes.

Genome-Wide Association Studies (GWAS) Data: Scans of genomes to find genetic markers associated with specific traits or diseases.

Single Nucleotide Polymorphism (SNP) Data: Study of genetic variations among individuals.

Copy Number Variation (CNV) Data: Study of variations in the number of copies of a particular gene.

Interactomics Data: Data related to the study of interactomes, which are the whole set of molecular interactions in a particular cell. This includes protein-protein interactions, protein-DNA interactions, and other types of interactions within the cell. It's crucial for understanding the functional networks and pathways in biological systems and how alterations in these networks can lead to disease.



# File formats 

Genomic Data (WGS and WES)

FASTA (.fasta, .fa): Stores nucleotide or peptide sequences.
FASTQ (.fastq): Stores sequences and their quality scores.
VCF (Variant Call Format, .vcf): Used for describing gene variants.

Transcriptomic Data

FASTQ (.fastq): For raw sequence data.
GTF/GFF (Gene Transfer Format/Generic Feature Format, .gtf/.gff): For annotating genes, exons, etc.

Proteomic Data

MASS Spectrometry Data Formats (.mzML, .raw, .mgf): For mass spectrometry data.
FASTA (.fasta): For peptide or protein sequences.

Metabolomic Data

NetCDF (Network Common Data Form, .cdf): For representing array-oriented data.
mzXML/mzData (.mzXML, .mzData): For mass spectrometry data.

Epigenetic Data

BAM/SAM (Binary Alignment/Map, Sequence Alignment/Map, .bam/.sam): For storing sequence data.
BED (Browser Extensible Data, .bed): For genomic intervals.

Phylogenetic Data

Newick Format (.newick, .nwk): For representing tree data.
NEXUS (.nexus): For richly annotated genetic sequence information.

Structural Biology Data

PDB (Protein Data Bank, .pdb): For 3D structures of large biological molecules.
mmCIF (macromolecular Crystallographic Information File, .cif): For crystallographic structures.

Microbiome Data

FASTA (.fasta): For sequence data.
BIOM (Biological Observation Matrix, .biom): For representing biological samples.

Population Genetics Data

PLINK (.ped, .map): For large-scale genotyping data.
VCF (.vcf): For variant information.

Clinical Data

HL7 (Health Level-7, various extensions): For health care and clinical data.
DICOM (Digital Imaging and Communications in Medicine, various extensions): For medical imaging.

Chemoinformatics Data

SDF (Structure-Data File, .sdf): For molecular structures and associated data.
SMILES (Simplified Molecular Input Line Entry System, .smi): For describing the structure of chemical molecules.

Pathway Data

SBML (Systems Biology Markup Language, .sbml): For systems biology models.
BioPAX (.owl): For biological pathway data.

GWAS Data

PLINK (.bed, .bim, .fam): Common format for GWAS data.
VCF (.vcf): For variant information.

SNP Data

VCF (.vcf): For describing SNP variants.
BED/BIM/FAM (PLINK formats): For SNP genotype data.

CNV Data

BED (.bed): For genomic intervals and CNV regions.
VCF (.vcf): Also used for CNV data.

Interactomics Data

PSI-MI XML (Proteomics Standards Initiative Molecular Interaction XML, .xml): For molecular interaction data.
MITAB (Molecular Interaction TAB delimited format, .mitab, .txt): For the tab-delimited format of molecular interactions.

# GWAS (Genome-Wide Association Studies) Data

Definition:
Genome-wide association Studies (GWAS) are research approaches used to identify genetic variants associated with specific traits, such as diseases, by scanning the genomes of many individuals. GWAS compares the DNA of participants with a trait or disease against those without to find SNPs (single nucleotide polymorphisms) that occur more frequently in those with the trait.

Example:
Imagine a GWAS investigating the genetic basis of diabetes. The study might identify an SNP located on chromosome 6 that is significantly associated with an increased risk of developing type 2 diabetes. This SNP could be represented on a Manhattan plot in the dashboard, with its position on the x-axis corresponding to its location on chromosome 6 and its -log10(p-value) on the y-axis indicating its level of significance.

# Sequencing Data
Definition:
Sequencing data refers to the comprehensive information obtained from sequencing the DNA or RNA of organisms. This data encompasses the precise order of nucleotides (adenine, thymine, cytosine, and guanine) in a genome or a specific gene, allowing for the detailed study of genetic variations, expression patterns, and more.

Example:
An example of sequencing data could be from a study examining the expression levels of genes across different cancer cell lines. The data might show that a particular gene, say Gene X, is highly expressed in one type of cancer cell line compared to normal tissue. This could be visualized on the dashboard through a heatmap, with rows representing different genes, columns representing different cell lines, and color intensity indicating the level of gene expression.








