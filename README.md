![ExploreASL](Design/ExploreASL_logoHeader.png)

<!-- ALL-CONTRIBUTORS-BADGE:START - Do not remove or modify this section -->
[![All Contributors](https://img.shields.io/badge/all_contributors-10-orange.svg?style=flat-square)](#contributors-) [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3905262.svg)](https://doi.org/10.5281/zenodo.3905262) [![View ExploreASL on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://nl.mathworks.com/matlabcentral/fileexchange/83203-exploreasl) ![GitHub tag (latest by date)](https://img.shields.io/github/v/tag/ExploreASL/ExploreASL) [![GitHub last commit](https://img.shields.io/github/last-commit/ExploreASL/Documentation?label=mkdocs)](https://exploreasl.github.io/Documentation/)
<!-- ALL-CONTRIBUTORS-BADGE:END -->

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

![ExploreASL Workflow](Design/WorkflowUpdate.png "Workflow of ExploreASL")

## Documentation

Additional information about ExploreASL can be found in the [Neuroimage paper]([https://pubmed.ncbi.nlm.nih.gov/32526385/) and on the [ExploreASL website](https://www.ExploreASL.org), including the walkthrough document and how-to videos. Further documentation is work in progress. For any help please contact the lead authors/developers at h.j.mutsaerts@amsterdamumc.nl or j.petr@hzdr.de.

## ExploreASL team

* [Henk Mutsaerts](mailto:h.j.mutsaerts@amsterdamumc.nl?subject=[GitHub]%20ExploreASL) - Co-creator, [Researcher](https://www.researchgate.net/profile/Henri-Mutsaerts)
* [Jan Petr](mailto:j.petr@hzdr.de?subject=[GitHub]%20ExploreASL) - Co-creator, [Researcher](https://www.researchgate.net/profile/Jan-Petr-2)
* [Michael Stritt](mailto:m.stritt@mediri.com?subject=[GitHub]%20ExploreASL) - Software developer, Research Associate, [ASPIRE](http://aspire-mri.eu/)
* [Paul Groot](mailto:p.f.c.groot@amsterdamumc.nl?subject=[GitHub]%20ExploreASL) - Software developer, IT specialist, [Researcher](https://www.researchgate.net/profile/Paul-Groot)
* [Pieter Vandemaele](mailto:pieter.vandemaele@gmail.com?subject=[GitHub]%20ExploreASL) - Developer Matlab [BIDS app](https://github.com/bids-standard)
* [Luigi Lorenzini](mailto:l.lorenzini@amsterdamumc.nl?subject=[GitHub]%20ExploreASL) - Developer [ExploreQC](https://www.hzdr.de/publications/Publ-31929)
* [Maurice Pasternak](mailto:maurice.pasternak@mail.utoronto.ca?subject=[GitHub]%20ExploreASL) - Developer [ExploreASL GUI](https://github.com/MauricePasternak/ExploreASL_GUI)
* [Mathijs Dijsselhof](mailto:m.b.dijsselhof@amsterdamumc.nl?subject=[GitHub]%20ExploreASL) - PhD student, [Cerebrovascular Age](https://sites.google.com/view/exploreasl/projects)
* [Beatriz Padrela](mailto:b.estevespadrela@amsterdamumc.nl?subject=[GitHub]%20ExploreASL) - PhD student, [BBB-ASL](https://sites.google.com/view/exploreasl/projects)
* [Sandeep Ganji](mailto:Sandeep.g.bio@gmail.com?subject=[GitHub]%20ExploreASL) - Developer integration Philips ISD, [Researcher](https://www.researchgate.net/profile/Sandeep-Ganji-3)
* [Patricia Clement](mailto:Patricia.Clement@ugent.be?subject=[GitHub]%20ExploreASL) - Developer ASL-BIDS & organizer, [Researcher](https://www.researchgate.net/profile/Patricia-Clement)

## Acknowledgments
This project is supported by the Dutch Heart Foundation (2020T049), the Eurostars-2 joint programme with co-funding from the European Union Horizon 2020 research and innovation programme ([ASPIRE E!113701](http://aspire-mri.eu/)), including the Netherlands Enterprise Agency (RvO), and by the EU Joint Program for Neurodegenerative Disease Research, including the Netherlands Organisation for health Research and Development and Alzheimer Nederland ([DEBBIE JPND2020-568-106](https://www.neurodegenerationresearch.eu/wp-content/uploads/2021/03/Project-DEBBIE.pdf)).

This project has previously received support from the following EU/EFPIA Innovative Medicines Initiatives (1 and 2) Joint Undertakings: [EPAD](http://ep-ad.org/) grant no. 115736, [AMYPAD](https://amypad.eu/) grant no. 115952 and [Amsterdam Neuroscience](https://www.amsterdamresearch.org/web/neuroscience/home.htm). The authors wish to thank the [COST-AID](https://asl-network.org/) (European Cooperation in Science and Technology - Arterial spin labeling Initiative in Dementia) Action BM1103 and the Open Source Initiative for Perfusion Imaging [(OSIPI)](https://www.osipi.org/) and the [ISMRM Perfusion Study groups](https://www.ismrm.org/study-groups/perfusion-mr/) for facilitating meetings for researchers to discuss the implementation of ExploreASL. The authors acknowledge Guillaume Flandin, Robert Dahnke, and Paul Schmidt for reviewing the structural module for its implementation of [SPM12](https://www.fil.ion.ucl.ac.uk/spm/software/spm12/), [CAT12](http://www.neuro.uni-jena.de/cat/), and [LST](https://www.applied-statistics.de/lst.html), respectively; Krzysztof Gorgolewksi for his advice on the [BIDS](https://bids.neuroimaging.io/) implementation; Jens Maus for help with MEX compilation; Cyril Pernet for providing the [SPM Univariate Plus](https://osf.io/wn3h8/) QC scripts.

## Contributors ✨

Thanks goes to these wonderful people ([emoji key](https://allcontributors.org/docs/en/emoji-key)):

<!-- ALL-CONTRIBUTORS-LIST:START - Do not remove or modify this section -->
<!-- prettier-ignore-start -->
<!-- markdownlint-disable -->
<table>
  <tr>
    <td align="center"><a href="http://www.ExploreASL.org"><img src="https://avatars0.githubusercontent.com/u/27774254?v=4" width="120px;" alt=""/><br /><sub><b>Henk Mutsaerts</b></sub></a><br /><a href="#creator-HenkMutsaerts" title="Mentor and Creator">👨‍🔬</a> <a href="#content-HenkMutsaerts" title="Content">🖋</a> <a href="https://github.com/ExploreASL/ExploreASL/commits?author=HenkMutsaerts" title="Code">💻</a></td>
    <td align="center"><a href="https://github.com/jan-petr"><img src="https://avatars0.githubusercontent.com/u/29886537?v=4" width="120px;" alt=""/><br /><sub><b>Jan Petr</b></sub></a><br /><a href="#creator-jan-petr" title="Mentor and Creator">👨‍🔬</a> <a href="#content-jan-petr" title="Content">🖋</a> <a href="https://github.com/ExploreASL/ExploreASL/commits?author=jan-petr" title="Code">💻</a></td>
    <td align="center"><a href="https://github.com/MichaelStritt"><img src="https://avatars0.githubusercontent.com/u/46593074?v=4" width="120px;" alt=""/><br /><sub><b>Michael Stritt</b></sub></a><br /><a href="https://github.com/ExploreASL/ExploreASL/commits?author=MichaelStritt" title="Code">💻</a> <a href="#content-MichaelStritt" title="Content">🖋</a> <a href="https://github.com/ExploreASL/ExploreASL/commits?author=MichaelStritt" title="Documentation">📖</a></td>
    <td align="center"><a href="https://github.com/MDijsselhof"><img src="https://avatars0.githubusercontent.com/u/75380250?v=4" width="120px;" alt=""/><br /><sub><b>Mathijs Dijsselhof</b></sub></a><br /><a href="#content-MDijsselhof" title="Content">🖋</a> <a href="#data-MDijsselhof" title="Data Acquisition & Management">🧠</a></td>
    <td align="center"><a href="https://github.com/BeatrizPadrela"><img src="https://avatars0.githubusercontent.com/u/73699072?v=4" width="120px;" alt=""/><br /><sub><b>Beatriz Padrela</b></sub></a><br /><a href="#content-BeatrizPadrela" title="Content">🖋</a> <a href="#data-BeatrizPadrela" title="Data Acquisition & Management">🧠</a></td>
  </tr>
  <tr>
    <td align="center"><a href="http://www.amsterdamumc.nl"><img src="https://avatars0.githubusercontent.com/u/18597189?v=4" width="120px;" alt=""/><br /><sub><b>Paul Groot</b></sub></a><br /><a href="https://github.com/ExploreASL/ExploreASL/commits?author=pfcgroot" title="Code">💻</a> <a href="#content-pfcgroot" title="Content">🖋</a></td>
    <td align="center"><a href="https://github.com/pvdemael"><img src="https://avatars1.githubusercontent.com/u/37624277?v=4" width="120px;" alt=""/><br /><sub><b>Pieter Vandemaele</b></sub></a><br /><a href="https://github.com/ExploreASL/ExploreASL/commits?author=pvdemael" title="Code">💻</a> <a href="#ideas-pvdemael" title="Ideas, Planning, & Feedback">🤔</a> <a href="#data-pvdemael" title="Data Acquisition & Management">🧠</a></td>
    <td align="center"><a href="https://github.com/MauricePasternak"><img src="https://avatars3.githubusercontent.com/u/57411571?v=4" width="120px;" alt=""/><br /><sub><b>MauricePasternak</b></sub></a><br /><a href="#gui-MauricePasternak" title="Graphical User Interface">📊</a> <a href="https://github.com/ExploreASL/ExploreASL/commits?author=MauricePasternak" title="Code">💻</a> <a href="#design-MauricePasternak" title="Design">🎨</a></td>
    <td align="center"><a href="https://github.com/patsycle"><img src="https://avatars0.githubusercontent.com/u/41481345?v=4" width="120px;" alt=""/><br /><sub><b>Patricia Clement</b></sub></a><br /> <a href="#data-patsycle" title="Data Acquisition & Management">🧠</a> <a href="#ideas-patsycle" title="Ideas, Planning, & Feedback">🤔</a> <a href="https://github.com/ExploreASL/ExploreASL/commits?author=patsycle" title="Documentation">📖</a> </td>
    <td align="center"><a href="https://github.com/sandeepganji"><img src="https://avatars0.githubusercontent.com/u/12124746?v=4" width="120px;" alt=""/><br /><sub><b>Sandeep Ganji</b></sub></a><br /><a href="#content-sandeepganji" title="Content">🖋</a> <a href="#ideas-sandeepganji" title="Ideas, Planning, & Feedback">🤔</a> <a href="#data-sandeepganji" title="Data Acquisition & Management">🧠</a></td>
  </tr>
  <tr>
    <td align="center"><a href="https://github.com/mcraig-ibme"><img src="https://avatars0.githubusercontent.com/u/26383586?v=4" width="120px;" alt=""/><br /><sub><b>Martin Craig</b></sub></a><br /><a href="#content-mcraig-ibme" title="Content">🖋</a> <a href="https://github.com/ExploreASL/ExploreASL/commits?author=mcraig-ibme" title="Code">💻</a> <a href="#data-mcraig-ibme" title="Data Acquisition & Management">🧠</a></td>
    <td align="center"><a href="https://github.com/DaveThoma5"><img src="https://avatars0.githubusercontent.com/u/3704113?v=4" width="120px;" alt=""/><br /><sub><b>DaveThoma5</b></sub></a><br /><a href="#ideas-DaveThoma5" title="Ideas, Planning, & Feedback">🤔</a> <a href="#data-DaveThoma5" title="Data Acquisition & Management">🧠</a></td>
    <td align="center"><a href="https://github.com/amahroo"><img src="https://avatars0.githubusercontent.com/u/48561261?v=4" width="120px;" alt=""/><br /><sub><b>Amnah Mahroo</b></sub></a><br /><a href="#ideas-amahroo" title="Ideas, Planning, & Feedback">🤔</a> <a href="#data-amahroo" title="Data Acquisition & Management">🧠</a></td>
    <td align="center"><a href="https://github.com/luislorenzini"><img src="https://avatars2.githubusercontent.com/u/57985241?v=4" width="120px;" alt=""/><br /><sub><b>luislorenzini</b></sub></a><br /><a href="https://github.com/ExploreASL/ExploreASL/commits?author=luislorenzini" title="Code">💻</a> <a href="#tool-luislorenzini" title="Tools">🔧</a></td>
    <td align="center"><a href="https://github.com/jozsait"><img src="https://avatars0.githubusercontent.com/u/19532128?v=4" width="120px;" alt=""/><br /><sub><b>jozsait</b></sub></a><br /><a href="https://github.com/ExploreASL/ExploreASL/commits?author=jozsait" title="Code">💻</a> <a href="#maintenance-jozsait" title="Maintenance">🚧</a></td>
  </tr>
</table>

<!-- markdownlint-enable -->
<!-- prettier-ignore-end -->
<!-- ALL-CONTRIBUTORS-LIST:END -->

This project follows the [all-contributors](https://github.com/all-contributors/all-contributors) specification. Contributions of any kind welcome!

## How to cite
The bare minimum of references is the [ExploreASL manuscript](https://www.sciencedirect.com/science/article/pii/S1053811920305176) and the used ExploreASL release, which you can find on Zenodo (e.g. [version 1.8.0](https://zenodo.org/record/5504194)).

The following provides an example as how to correctly cite ExploreASL and its third-party tools. The versions of the included third-party tools are described in [CHANGES.md](https://github.com/ExploreASL/ExploreASL/blob/main/CHANGES.md) for each ExploreASL release.

>The data were analysed using ExploreASL `ref1` version x.x.x `ref2`, including SPM12 version xxxx `ref3`, CAT12 version xxxx`ref4`, and LST version x.x.x`ref5`. This Matlab-based software was used with Matlab (MathWorks, MA, USA) version x.x (yearx)`ref6`.

### References

The release numbers of ExploreASL (e.g. `1.9.0`) follow [semantic versioning](https://semver.org/).

1. The [ExploreASL paper](https://www.sciencedirect.com/science/article/pii/S1053811920305176), describing the full pipeline and decisions for processing steps.
2. The Zenodo DOI for the actual ExploreASL release used to analyse the data (e.g. the [up-to-date release](https://doi.org/10.5281/zenodo.3905262)).
3. The SPM12 references [Ashburner, 2012](https://doi.org/10.1016/j.neuroimage.2011.10.025) & [Flandin and Friston, 2008](https://doi.org/10.4249/scholarpedia.6232). Note that the SPM version (e.g. `7219`) is adapted and extended for use with ExploreASL.
4. The CAT12 reference [Gaser, 2009](https://doi.org/10.1016/S1053-8119(09)71151-6). Note that the CAT12 version (e.g. `1364`) is adapted for use with ExploreASL.
5. The LST references [Schmidt, 2017](https://www.statistical-modelling.de/LST_documentation.pdf) & [de Sitter, 2017](https://doi.org/10.1016/j.neuroimage.2017.09.011). Note that the LST version (e.g. `2.0.15`) is adapted for use with ExploreASL.
6. Matlab publishes a [release](https://www.mathworks.com/help/matlab/release-notes.html) twice yearly. You can provide the release number (e.g. `9.4`) or year number (e.g. `2018a`), or both.
