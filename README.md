![ExploreASL](https://github.com/ExploreASL/ExploreASL/blob/develop/ExploreASL_logoSmall.png)

# ExploreASL

<!-- ALL-CONTRIBUTORS-BADGE:START - Do not remove or modify this section -->
[![All Contributors](https://img.shields.io/badge/all_contributors-10-orange.svg?style=flat-square)](#contributors-)
<!-- ALL-CONTRIBUTORS-BADGE:END -->

<!-- ZENODO-BADGE:START - Do not remove or modify this section -->
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3905262.svg)](https://doi.org/10.5281/zenodo.3905262)
<!-- ZENODO-BADGE:END -->

* [ExploreASL](#exploreasl)
	* [Description](#description)
	* [Installation](#installation)
	* [Workflow](#workflow)
	* [Documentation](#documentation)
	* [ExploreASL team](#exploreasl-team)
	* [Acknowledgments](#acknowledgments)
	* [How to cite](#how-to-cite)

## Description

**ExploreASL** is a pipeline and toolbox for image processing and statistics of arterial spin labeling perfusion **MR** images. It is designed as a **multi-OS**, open source, collaborative framework that facilitates cross-pollination between image processing method developers and clinical investigators.

The software provides a complete head-to-tail approach that runs fully automatically, encompassing all necessary tasks from data import and structural segmentation, registration and normalization, up to **CBF** quantification. In addition, the software package includes and quality control (**QC**) procedures and region-of-interest (**ROI**) as well as voxel-wise analysis on the extracted data. To-date, **ExploreASL** has been used for processing ~10000 **ASL** datasets from all major **MRI** vendors and **ASL** sequences, and a variety of patient populations, representing ~30 studies. The ultimate goal of **ExploreASL** is to combine data from multiple studies to identify disease related perfusion patterns that may prove crucial in using **ASL** as a diagnostic tool and enhance our understanding of the interplay of perfusion and structural changes in neurodegenerative pathophysiology. 

Additionally, this (semi-)automatic pipeline allows us to minimize manual intervention, which increases the reproducibility of studies. 

## Installation

To use **ExploreASL** within Matlab, you can download a stable release version from the GitHub releases section. Navigate within Matlab to the **ExporeASL** directory, to make **ExploreASL** the current working directory. To start ExploreASL from Matlab, type:

```
ExploreASL
```

## Workflow

![ExploreASL Workflow](https://www.researchgate.net/profile/Andrew_Robertson7/publication/337328693/figure/fig1/AS:826578854481921@1574083164220/Schematic-diagram-of-ExploreASL-processing-steps-Steps-marked-with-a-are-optional.ppm "Workflow of ExploreASL")

## Documentation

Additional information about ExploreASL can be found in the [Neuroimage paper]([https://pubmed.ncbi.nlm.nih.gov/32526385/) and on the [ExploreASL](www.ExploreASL.org) website, including the walkthrough document and how-to videos. Further documentation is work in progress. For any help please contact the lead authors/developers at h.j.mutsaerts@amsterdamumc.nl or j.petr@hzdr.de.

## ExploreASL team

* [Henk Mutsaerts](h.j.mutsaerts@amsterdamumc.nl) - co-creator
* [Jan Petr](j.petr@hzdr.de) - co-creator
* [Michael Stritt](stritt.michael@gmail.com) - ASPIRE PhD
* [Paul Groot](p.f.c.groot@amsterdamumc.nl) - developer backbone, IT specialist
* [Pieter Vandemaele](pieter.vandemaele@gmail.com) - developer Matlab BIDS app
* [Luigi Lorenzini](l.lorenzini@amsterdamumc.nl) - developer ExploreQC
* [Maurice Pasternak](maurice.pasternak@mail.utoronto.ca) - developer GUI
* [Sandeep Ganji](Sandeep.g.bio@gmail.com) - developer integration Philips ISD
* [Patricia Clement](Patricia.Clement@ugent.be) - developer ASL-BIDS & organizer

## Acknowledgments
This project has received support from the following EU/EFPIA Innovative Medicines Initiatives (1 and 2) Joint Undertakings: [EPAD](http://ep-ad.org/) grant no. 115736, [AMYPAD](https://amypad.eu/) grant no. 115952. Additionally, this work received support from the EU-EFPIA Innovative Medicines Initiatives Joint Undertaking (grant No 115952), and [Amsterdam Neuroscience](https://www.amsterdamresearch.org/web/neuroscience/home.htm). The authors wish to thank the [COST-AID](https://asl-network.org/) (European Cooperation in Science and Technology - Arterial spin labeling Initiative in Dementia) Action BM1103 and the Open Source Initiative for Perfusion Imaging [(OSIPI)](https://www.osipi.org/) and the [ISMRM Perfusion Study groups](https://www.ismrm.org/study-groups/perfusion-mr/) for facilitating meetings for researchers to discuss the implementation of ExploreASL. The authors acknowledge Guillaume Flandin, Robert Dahnke, and Paul Schmidt for reviewing the structural module for its implementation of [SPM12](https://www.fil.ion.ucl.ac.uk/spm/software/spm12/), [CAT12](http://www.neuro.uni-jena.de/cat/), and [LST](https://www.applied-statistics.de/lst.html), respectively; Krzysztof Gorgolewksi for his advice on the [BIDS](https://bids.neuroimaging.io/) implementation; Jens Maus for help with MEX compilation; Cyril Pernet for providing the [SPM Univariate Plus](https://osf.io/wn3h8/) QC scripts.

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
    <td align="center"><a href="http://www.amsterdamumc.nl"><img src="https://avatars0.githubusercontent.com/u/18597189?v=4" width="100px;" alt=""/><br /><sub><b>Paul Groot</b></sub></a><br /><a href="https://github.com/ExploreASL/ExploreASL/commits?author=pfcgroot" title="Code">💻</a> <a href="#content-pfcgroot" title="Content">🖋</a></td>
    <td align="center"><a href="https://github.com/pvdemael"><img src="https://avatars1.githubusercontent.com/u/37624277?v=4" width="100px;" alt=""/><br /><sub><b>Pieter Vandemaele</b></sub></a><br /><a href="https://github.com/ExploreASL/ExploreASL/commits?author=pvdemael" title="Code">💻</a> <a href="#ideas-pvdemael" title="Ideas, Planning, & Feedback">🤔</a></td>
    <td align="center"><a href="https://github.com/luislorenzini"><img src="https://avatars2.githubusercontent.com/u/57985241?v=4" width="100px;" alt=""/><br /><sub><b>luislorenzini</b></sub></a><br /><a href="https://github.com/ExploreASL/ExploreASL/commits?author=luislorenzini" title="Code">💻</a> <a href="#tool-luislorenzini" title="Tools">🔧</a></td>
    <td align="center"><a href="https://github.com/MauricePasternak"><img src="https://avatars3.githubusercontent.com/u/57411571?v=4" width="100px;" alt=""/><br /><sub><b>MauricePasternak</b></sub></a><br /><a href="https://github.com/ExploreASL/ExploreASL/commits?author=MauricePasternak" title="Code">💻</a> <a href="#design-MauricePasternak" title="Design">🎨</a></td>
  </tr>
  <tr>
    <td align="center"><a href="https://github.com/sandeepganji"><img src="https://avatars0.githubusercontent.com/u/12124746?v=4" width="100px;" alt=""/><br /><sub><b>Sandeep Ganji</b></sub></a><br /><a href="#content-sandeepganji" title="Content">🖋</a> <a href="#ideas-sandeepganji" title="Ideas, Planning, & Feedback">🤔</a></td>
    <td align="center"><a href="https://github.com/jozsait"><img src="https://avatars0.githubusercontent.com/u/19532128?v=4" width="100px;" alt=""/><br /><sub><b>jozsait</b></sub></a><br /><a href="https://github.com/ExploreASL/ExploreASL/commits?author=jozsait" title="Code">💻</a> <a href="#maintenance-jozsait" title="Maintenance">🚧</a></td>
    <td align="center"><a href="https://github.com/DaveThoma5"><img src="https://avatars0.githubusercontent.com/u/3704113?v=4" width="100px;" alt=""/><br /><sub><b>DaveThoma5</b></sub></a><br /><a href="#ideas-DaveThoma5" title="Ideas, Planning, & Feedback">🤔</a></td>
  </tr>
</table>

<!-- markdownlint-enable -->
<!-- prettier-ignore-end -->
<!-- ALL-CONTRIBUTORS-LIST:END -->

This project follows the [all-contributors](https://github.com/all-contributors/all-contributors) specification. Contributions of any kind welcome!

## How to cite
The bare minimum of references is the [ExploreASL manuscript](https://www.sciencedirect.com/science/article/pii/S1053811920305176) and the used ExploreASL release, which you can find on Zenodo (e.g. [version 1.3.0](https://zenodo.org/record/4095518#.X4rTL5P7Rts)).

The following provides an example as how to correctly cite ExploreASL and its third-party tools. The versions of the included third-party tools are described in [CHANGES.md](https://github.com/ExploreASL/ExploreASL/blob/master/CHANGES.md) for each ExploreASL release.

>The data were analysed using ExploreASL `ref1` version x.x.x `ref2`, including SPM12 version xxxx `ref3`, CAT12 version xxxx`ref4`, and LST version x.x.x`ref5`. This Matlab-based software was used with Matlab (MathWorks, MA, USA) version x.x (yearx)`ref6`.

* Ref1: the ExploreASL paper, describing the full pipeline and decisions for processing steps: https://www.sciencedirect.com/science/article/pii/S1053811920305176
* Ref2: the Zenodo DOI for the actual ExploreASL release used to analyse the data. The release numbers (e.g. 1.3.0) follow [semantic versioning](https://semver.org/).
* Ref3: SPM12 references: https://www.sciencedirect.com/science/article/pii/S1053811920305176#bib14 & https://www.sciencedirect.com/science/article/pii/S1053811920305176#bib53. Note that the SPM version (e.g. 7219) is adapted and extended for use with ExploreASL.
* Ref4: CAT12 reference: https://www.sciencedirect.com/science/article/pii/S1053811920305176#bib55. Note that the CAT12 version (e.g. 1364) is adapted for use with ExploreASL.
* Ref5: LST reference: https://www.sciencedirect.com/science/article/pii/S1053811920305176#bib118. Note that the LST version (e.g. 2.0.15) is adapted for use with ExploreASL.
* Ref6: Matlab publishes a release twice yearly, which can be reviewed here: https://www.mathworks.com/products/compiler/matlab-runtime.html. You can provide the release number (e.g. 9.4) or year number (e.g. 2018a), or both.
