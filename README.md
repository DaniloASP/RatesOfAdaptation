# Estimating the rates of adaptive evolution in fungal species using whole-genome sequences

[![MIT License][license-shield]][license-url]
[![LinkedIn][linkedin-shield]][linkedin-url]


# Content

-    [Description](#description)
-    [Installation](#installation)
-    [Usage](#usage)

# Description
These scripts are used for the calculation of the rates of adaptation from raw Illumina reads. The rates of adaptation are estimated by means of:
* $\alpha$ (alpha): Proportion of non-synonymous mutations fixed by positive selection 
* $\omega$ A (omegaA): Rate of adaptive non-synonymous substitutions
* $\omega$ NA (omegaNA): Rate of non-adaptive non-synonymous substitutions

# Installation
These scripts use a series of tools that you might need depending on the step you want to run. Some important tools are:

**mafft** is used to add the outgroup sequence to the aligned fasta file and perform a new alignment. Instruction can be found in https://mafft.cbrc.jp/alignment/software/

**MACSE** is used to polish the alignment output of maft. Installation instructions can be found in https://github.com/ranwez/MACSE_V2_PIPELINES

**Maffilter** is used for the extraction of CDS/exon regions from VCF files. For details see https://jydu.github.io/maffilter/

**BppPopstat** is used for the calculation of kappa (the transition over transversion ratio) and quantification of polymorphism and divergence mutations across gene alignments. Installation instruction can be found in https://github.com/BioPP/bppsuite or [here](https://github.com/DaniloASP/RatesOfAdaptation/blob/main/4_bppPopStat/Installation_bppPopStat.sh).

**CNVnator** is used to identify deleted regions across the genome. Installation instructions can be found in https://github.com/abyzovlab/CNVnator

**LDhat** is used to calculate the population recombination rate ($\rho$). Installation and look-up tables can be found in https://github.com/auton1/LDhat

**ANGSD** is used to estimate Watterson $\theta$ for the selection of the most appropriate look-up table from LDhat. Installation instructions can be found in https://github.com/ANGSD/angsd

**GRAPES** is used to estimate the rates of adaptive evolution. Installation instructions can be found in https://github.com/BioPP/grapes

# Usage

**Data preparation:** Raw illumina reads are trimmed, filtered for quality and mapped onto a reference genome. After hardfiltering and sanity checks, a consensus genome is extracted from each sample in the VCF, and Maffilter is used to extract the coding region of reach consensus genome. A custom python script is then used to merge all conding regions per gene. Scripts in [1_DataPreparation](https://github.com/DaniloASP/RatesOfAdaptation/tree/main/1_DataPreparation).

**Gene alignment:** Orthologous genes will be identified among the reference isolates. For the population data, individual genes files will extracted into fasta files containing isolate's sequence and the outgroup sequence will be added. Mafft and MACSE will be used for the alignment and polishment, and custom python scripts will be used for quality control. Scripts in [2_GeneAlignment](https://github.com/DaniloASP/RatesOfAdaptation/tree/main/2_GeneAlignment).

**Protein functional analysis**: Using interpro scan and effectorP to classify proteins into non-secreted, secreted but not effector and predicted effectors. Scripts in [3_Protein_FunctionalAnalysis](https://github.com/DaniloASP/RatesOfAdaptation/tree/main/3_Protein_FunctionalAnalysis).

**bppPopStat:** Calculation of kappa, polymorphism and substitutions counts are performed using this tool on each gene individually. A first run is performed for the estimation of the median kappa across all genes. Here, the distribution of kappa must be unimodal (see Fig1). Next, kappa is fixed to the previously determined median and the counts of mutations are performed. [4_bppPopStat](https://github.com/DaniloASP/RatesOfAdaptation/tree/main/4_bppPopStat).

<img src="https://github.com/DaniloASP/RatesOfAdaptation/blob/main/4_bppPopStat/Zymoseptoriatritici.485.4935genes.kappa3.05.core.jpeg" alt="drawing" width="600"/>

Fig1. Example of kappa estimations in the fungus *Zymoseptoria tritici*


<!-- MARKDOWN LINKS & IMAGES -->
[license-shield]: https://img.shields.io/github/license/othneildrew/Best-README-Template.svg?style=for-the-badge
[license-url]: https://github.com/DaniloASP/RatesOfAdaptation/blob/main/LICENSE
[linkedin-shield]: https://img.shields.io/badge/-LinkedIn-black.svg?style=for-the-badge&logo=linkedin&colorB=555
[linkedin-url]: https://linkedin.com/in/das-pereira

