![ExploreASL](Design/ExploreASL_logoHeader.png)

<!-- ALL-CONTRIBUTORS-BADGE:START - Do not remove or modify this section -->
[![All Contributors](https://img.shields.io/badge/all_contributors-10-orange.svg?style=flat-square)](#contributors-) [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3905262.svg)](https://doi.org/10.5281/zenodo.3905262) [![View ExploreASL on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://nl.mathworks.com/matlabcentral/fileexchange/83203-exploreasl) ![GitHub tag (latest by date)](https://img.shields.io/github/v/tag/ExploreASL/ExploreASL) [![GitHub last commit](https://img.shields.io/github/last-commit/ExploreASL/Documentation?label=mkdocs)](https://exploreasl.github.io/Documentation/)
<!-- ALL-CONTRIBUTORS-BADGE:END -->

## Description

**ExploreASL** is a pipeline and toolbox for image processing and statistics of arterial spin labeling perfusion **MR** images. It is designed as a **multi-OS**, open-source, collaborative framework that facilitates cross-pollination between image processing method developers and clinical investigators.

The software provides a complete head-to-tail approach that runs fully automatically, encompassing all necessary tasks from data import and structural segmentation, registration, and normalization, up to **CBF** quantification. In addition, the software package includes quality control (**QC**) procedures and region-of-interest (**ROI**) as well as voxel-wise analysis of the extracted data. To date, **ExploreASL** has been used for processing ~10000 **ASL** datasets from all major **MRI** vendors and **ASL** sequences and a variety of patient populations, representing ~30 studies. The ultimate goal of **ExploreASL** is to combine data from multiple studies to identify disease-related perfusion patterns that may prove crucial in using **ASL** as a diagnostic tool and enhance our understanding of the interplay of perfusion and structural changes in neurodegenerative pathophysiology. 

Additionally, this (semi-)automatic pipeline allows us to minimize manual intervention, which increases the reproducibility of studies. 

## Documentation

Reference manual and tutorials for each ExploreASL version are found on the [GitHub website](https://exploreasl.github.io/Documentation). A general description of ExploreASL is in the [Neuroimage paper]([https://pubmed.ncbi.nlm.nih.gov/32526385/). Additional resources are on the [ExploreASL website](https://www.ExploreASL.org) including the walkthrough document and how-to videos, but these are not regularly updated with new versions. For any help please contact the lead authors/developers at h.j.mutsaerts@amsterdamumc.nl or j.petr@hzdr.de.

## Installation

To use **ExploreASL** within Matlab, you can download a stable release version from the [GitHub releases section](https://github.com/ExploreASL/ExploreASL/releases) or from [Zenodo](https://zenodo.org/record/3905263). Navigate within Matlab to the **ExporeASL** directory, to make **ExploreASL** the current working directory. To start ExploreASL from Matlab, type:

```
ExploreASL
```

## Workflow

![ExploreASL Workflow](Design/WorkflowUpdate.png "Workflow of ExploreASL")

## Acknowledgments
This project is supported by the Dutch Heart Foundation (2020T049), the Eurostars-2 joint programme with co-funding from the European Union Horizon 2020 research and innovation programme ([ASPIRE E!113701](http://aspire-mri.eu/)), including the Netherlands Enterprise Agency (RvO), and by the EU Joint Program for Neurodegenerative Disease Research, including the Netherlands Organisation for health Research and Development and Alzheimer Nederland ([DEBBIE JPND2020-568-106](https://www.neurodegenerationresearch.eu/wp-content/uploads/2021/03/Project-DEBBIE.pdf)).

This project has previously received support from the following EU/EFPIA Innovative Medicines Initiatives (1 and 2) Joint Undertakings: [EPAD](http://ep-ad.org/) grant no. 115736, [AMYPAD](https://amypad.eu/) grant no. 115952 and [Amsterdam Neuroscience](https://www.amsterdamresearch.org/web/neuroscience/home.htm). The authors wish to thank the [COST-AID](https://asl-network.org/) (European Cooperation in Science and Technology - Arterial spin labeling Initiative in Dementia) Action BM1103 and the Open Source Initiative for Perfusion Imaging [(OSIPI)](https://www.osipi.org/) and the [ISMRM Perfusion Study groups](https://www.ismrm.org/study-groups/perfusion-mr/) for facilitating meetings for researchers to discuss the implementation of ExploreASL. The authors acknowledge Guillaume Flandin, Robert Dahnke, and Paul Schmidt for reviewing the structural module for its implementation of [SPM12](https://www.fil.ion.ucl.ac.uk/spm/software/spm12/), [CAT12](http://www.neuro.uni-jena.de/cat/), and [LST](https://www.applied-statistics.de/lst.html), respectively; Krzysztof Gorgolewksi for his advice on the [BIDS](https://bids.neuroimaging.io/) implementation; Jens Maus for help with MEX compilation; Cyril Pernet for providing the [SPM Univariate Plus](https://osf.io/wn3h8/) QC scripts.

## How to cite
The following provides an example as how to correctly cite ExploreASL and its third-party tools. The versions of the included third-party tools are described in [CHANGES.md](https://github.com/ExploreASL/ExploreASL/blob/main/CHANGES.md) for each ExploreASL release. The bare minimum of references (`refs`) are `ref1` and `ref2`.

>The data were analysed using ExploreASL `ref1` version x.x.x `ref2`, including SPM12 version xxxx `ref3`, CAT12 version xxxx`ref4`, and LST version x.x.x`ref5`. This Matlab-based software was used with Matlab (MathWorks, MA, USA) version x.x (yearx)`ref6`.

### References

The release numbers of ExploreASL (e.g. `1.9.0`) follow [semantic versioning](https://semver.org/).

1. The [ExploreASL paper](https://www.sciencedirect.com/science/article/pii/S1053811920305176), describing the full pipeline and decisions for processing steps.
2. The Zenodo DOI for the actual ExploreASL release used to analyse the data, e.g. the [latest release](https://doi.org/10.5281/zenodo.3905262)).
3. The SPM12 references [Ashburner, 2012](https://doi.org/10.1016/j.neuroimage.2011.10.025) & [Flandin and Friston, 2008](https://doi.org/10.4249/scholarpedia.6232). Note that the SPM version (e.g. `7219`) is adapted and extended for use with ExploreASL.
4. The CAT12 reference [Gaser, 2009](https://doi.org/10.1016/S1053-8119(09)71151-6). Note that the CAT12 version (e.g. `1364`) is adapted for use with ExploreASL.
5. The LST references [Schmidt, 2017](https://www.statistical-modelling.de/LST_documentation.pdf) & [de Sitter, 2017](https://doi.org/10.1016/j.neuroimage.2017.09.011). Note that the LST version (e.g. `2.0.15`) is adapted for use with ExploreASL.
6. Matlab publishes a [release](https://www.mathworks.com/help/matlab/release-notes.html) twice yearly. You can provide the release number (e.g. `9.4`) or year number (e.g. `2018a`), or both.
