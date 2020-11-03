# ExploreASL
<!-- ALL-CONTRIBUTORS-BADGE:START - Do not remove or modify this section -->
[![All Contributors](https://img.shields.io/badge/all_contributors-4-orange.svg?style=flat-square)](#contributors-)
<!-- ALL-CONTRIBUTORS-BADGE:END -->

## Description

**ExploreASL** is a pipeline and toolbox for image processing and statistics of arterial spin labeling perfusion **MR** images. It is designed as a **multi-OS**, open source, collaborative framework that facilitates cross-pollination between image processing method developers and clinical investigators.

The software provides a complete head-to-tail approach that runs fully automatically, encompassing all necessary tasks from data import and structural segmentation, registration and normalization, up to **CBF** quantification. In addition, the software package includes and quality control (**QC**) procedures and region-of-interest (**ROI**) as well as voxel-wise analysis on the extracted data. To-date, **ExploreASL** has been used for processing ~10000 **ASL** datasets from all major **MRI** vendors and **ASL** sequences, and a variety of patient populations, representing ~30 studies. The ultimate goal of **ExploreASL** is to combine data from multiple studies to identify disease related perfusion patterns that may prove crucial in using **ASL** as a diagnostic tool and enhance our understanding of the interplay of perfusion and structural changes in neurodegenerative pathophysiology. 

Additionally, this (semi-)automatic pipeline allows us to minimize manual intervention, which increases the reproducibility of studies. 

## Highlighted Features

## Features

## Installation

To use **ExploreASL** within Matlab, you can download a stable release version from the GitHub releases section. Navigate within Matlab to the **ExporeASL** directory, to make **ExploreASL** the current working directory. To start ExploreASL from Matlab, type:

```
ExploreASL
```


## Workflow

![ExploreASL Workflow](https://www.researchgate.net/profile/Andrew_Robertson7/publication/337328693/figure/fig1/AS:826578854481921@1574083164220/Schematic-diagram-of-ExploreASL-processing-steps-Steps-marked-with-a-are-optional.ppm "Workflow of ExploreASL")

## Documentation
Additional information about ExploreASL can be found in the [Neuroimage paper]([https://pubmed.ncbi.nlm.nih.gov/32526385/) and on the [ExploreASL](www.ExploreASL.org) website, including the walkthrough document and how-to videos. Further documentation is work in progress. For any help please contact the lead authors/developers at h.j.mutsaerts@amsterdamumc.nl or j.petr@hzdr.de.

## Creators

Please contact the co-creators for more information:

* Henk-Jan Mutsaerts *HenkJanMutsaerts@Gmail.com*
* Jan Petr *j.petr@hzdr.de*

## Developer team
* Michael Stritt *stritt.michael@gmail.com*
* Paul Groot *p.f.c.groot@amsterdamumc.nl*
* Pieter Vandemaele *pieter.vandemaele@gmail.com*
* Luigi Lorenzini *l.lorenzini@amsterdamumc.nl*
* Maurice Pasternak *maurice.pasternak@mail.utoronto.ca*
* Sandeep Ganji *Sandeep.g.bio@gmail.com*

## Acknowledgement
This project has received support from the following EU/EFPIA Innovative Medicines Initiatives (1 and 2) Joint Undertakings: [EPAD](http://ep-ad.org/) grant no. 115736, [AMYPAD](https://amypad.eu/) grant no. 115952. Additionally, this work received support from the EU-EFPIA Innovative Medicines Initiatives Joint Undertaking (grant No 115952), and [Amsterdam Neuroscience](https://www.amsterdamresearch.org/web/neuroscience/home.htm). The authors wish to thank the [COST-AID](https://asl-network.org/) (European Cooperation in Science and Technology - Arterial spin labeling Initiative in Dementia) Action BM1103 and the Open Source Initiative for Perfusion Imaging [(OSIPI)](https://www.osipi.org/) and the [ISMRM Perfusion Study groups](https://www.ismrm.org/study-groups/perfusion-mr/)for facilitating meetings for researchers to discuss the implementation of ExploreASL. The authors acknowledge Guillaume Flandin, Robert Dahnke, and Paul Schmidt for reviewing the structural module for its implementation of [SPM12](https://www.fil.ion.ucl.ac.uk/spm/software/spm12/), [CAT12](http://www.neuro.uni-jena.de/cat/), and [LST](https://www.applied-statistics.de/lst.html), respectively; Krzysztof Gorgolewksi for his advice on the [BIDS](https://bids.neuroimaging.io/) implementation; Jens Maus for help with MEX compilation; Cyril Pernet for providing the [SPM Univariate Plus](https://osf.io/wn3h8/) QC scripts.

## Contributors ✨

Thanks goes to these wonderful people ([emoji key](https://allcontributors.org/docs/en/emoji-key)):

<!-- ALL-CONTRIBUTORS-LIST:START - Do not remove or modify this section -->
<!-- prettier-ignore-start -->
<!-- markdownlint-disable -->
<table>
  <tr>
    <td align="center"><a href="http://www.ExploreASL.org"><img src="https://avatars0.githubusercontent.com/u/27774254?v=4" width="100px;" alt=""/><br /><sub><b>Henk Mutsaerts</b></sub></a><br /><a href="#design-HenkMutsaerts" title="Design">🎨</a> <a href="#content-HenkMutsaerts" title="Content">🖋</a></td>
    <td align="center"><a href="https://github.com/jan-petr"><img src="https://avatars0.githubusercontent.com/u/29886537?v=4" width="100px;" alt=""/><br /><sub><b>Jan Petr</b></sub></a><br /><a href="#content-jan-petr" title="Content">🖋</a> <a href="#design-jan-petr" title="Design">🎨</a></td>
    <td align="center"><a href="https://github.com/MichaelStritt"><img src="https://avatars0.githubusercontent.com/u/46593074?v=4" width="100px;" alt=""/><br /><sub><b>Michael Stritt</b></sub></a><br /><a href="#content-MichaelStritt" title="Content">🖋</a> <a href="https://github.com/ExploreASL/ExploreASL/commits?author=MichaelStritt" title="Documentation">📖</a></td>
    <td align="center"><a href="https://github.com/sandeepganji"><img src="https://avatars0.githubusercontent.com/u/12124746?v=4" width="100px;" alt=""/><br /><sub><b>Sandeep Ganji</b></sub></a><br /><a href="#content-sandeepganji" title="Content">🖋</a> <a href="#ideas-sandeepganji" title="Ideas, Planning, & Feedback">🤔</a></td>
  </tr>
</table>

<!-- markdownlint-enable -->
<!-- prettier-ignore-end -->
<!-- ALL-CONTRIBUTORS-LIST:END -->

This project follows the [all-contributors](https://github.com/all-contributors/all-contributors) specification. Contributions of any kind welcome!