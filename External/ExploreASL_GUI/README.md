# Explore ASL GUI

> Complementary GUI to assist Arterial Spin Labelling analysis by ExploreASL

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
![Version](https://img.shields.io/badge/Version-0.2.0-yellow)
![PythonVersions](https://img.shields.io/badge/Python-3.7%20%7C%203.8-green)

## Graphics and Examples

> Drag and Drop your variables of interest, then view the subject to which a datapoint corresponds!

![MRI_Viewer_movie](https://github.com/MauricePasternak/ExploreASL_GUI/blob/master/github_media/MRI_viewer.gif)

---
## Table of Contents
- [Installation](#installation)
- [Features](#features)
- [Team](#exploreasl-team)
- [Support](#support)
- [FAQ](#faq)
- [License](#license)

---
## Installation

To be updated soon

---
## Features

The graphical user interface is split into 4 main modules existing as separate windows linked together by a central hub from which they can be launched:

![Image of MainWindow](https://github.com/MauricePasternak/ExploreASL_GUI/blob/master/github_media/MainWindow.png)

The first module most users will utilize is the Importer Module. The user defines the source directory (root) of their dataset and which eventually contains dicom files at terminal directories within the folder structure. This module then converts the dicoms therein into nifti volumes within a BIDS-compliant data structure, [utilizing Chris Rorden's dcm2niix executable in the process](https://github.com/rordenlab/dcm2niix):

![Image of Importer](https://github.com/MauricePasternak/ExploreASL_GUI/blob/master/github_media/Importer.png)

The main ExploreASL program has over 40 adjustable parameters to account for ASL acquisition variability. Defining these parameters is simplified through the ParmsMaker Module and the resulting DataPar.json file is automatically saved to the indicated study directory. In addition, existing DataPar.json files can be re-loaded, edited, and then re-saved back into the study directory.

![Image of ParmsMaker](https://github.com/MauricePasternak/ExploreASL_GUI/blob/master/github_media/ParmsMaker.png)

Once data is imported and parameters defined, the user will be prepared to run the main program. The Executor Module allows for multi-processing analysis of multiple studies and/or multiple subjects per study simultaneously.

![Image of Executor](https://github.com/MauricePasternak/ExploreASL_GUI/blob/master/github_media/Executor.png)

Finally, with the data analyzed, the user may feel free to visualize the processed dataset. The Post Processing & Visualization Module allows for users to interactively drag & drop all loaded data to generate publication-quality plots and load in both functional & structural MRI scans based on user clicks on specific datapoints.

![Image of PostProc](https://github.com/MauricePasternak/ExploreASL_GUI/blob/master/github_media/PostProc.png)

From start to finish, this GUI is designed to streamline the processing and eventual publication of arterial spin labelling image data.

This program is intended solely for research purposes and **should not** be used in any way, shape, or form for clinical application or in a clinical setting.

---
## Creator

For more information on the GUI, troubleshooting, and suggestions for improvements, please don't hesitate to contact:
- Maurice Pasternak; maurice.pasternak@utoronto.ca

---
## ExploreASL Team

For questions or concerns with the underlying ExploreASL program, the following individuals compromise the development team and may be contacted:

- Henk-Jan Mutsaerts; HenkJanMutsaerts@Gmail.com (ExploreASL creator)
- Jan Petr; j.petr@hzdr.de (ExploreASL creator)
- Michael Stritt; stritt.michael@gmail.com
- Paul Groot; p.f.c.groot@amsterdamumc.nl
- Pieter Vandemaele; pieter.vandemaele@gmail.com
- Maurice Pasternak; maurice.pasternak@utoronto.ca
- Luigi Lorenzini; l.lorenzini@amsterdamumc.nl
- Sandeep Ganji; Sandeep.g.bio@gmail.com

---
## Support

For GUI-related questions, please don't hesitate to drop an email at: maurice.pasternak@utoronto.ca

For more information on the main program, click on the following: [CLICK ME!](https://sites.google.com/view/exploreasl)

For the Github page of main program that this GUI interfaces with, click on the following: [CLICK ME!](https://github.com/ExploreASL/ExploreASL)

---
## FAQ

> **Q: What operating systems is your program compatible with?**

A: The underlying toolkit utilized in the creation of this program is PySide2, a Python wrapper around the cross-platform C++ software Qt. As such, the program should be compatible with the most common operating systems, including:
- Windows 7/Vista/10
- Mac OS X
- Linux

> **Q: Do I need MATLAB in order to run this program?**

A: Unfortunately, as ExploreASL itself has not been compiled into a binary executable, the GUI currently remains reliant on the user having an installed & activated version of MATLAB on their machine. Upon ExploreASL becoming a separate executable, this GUI will no longer require MATLAB to be present.

> **Q: The number of cores listed in your Executor module is off. My machine has ___ cores but you list ___ instead.**

A: This is primarily an issue on older and/or lower-binned CPUs that do not have hyper-threading (Intel) or simultaneous multi-threading (AMD) enabled, as the GUI follows the old adage "better safe than sorry". The number of avaliable logical cores avaliable is divided by a factor of 2. For CPUs with hyperthreading, this is a non-issue, as the resulting number will be equal to the number of physical cores. This is acceptable because the foundation behind ExploreASL - SPM12, is itself capable of utilizing both logical cores within each assigned physical core. For machines within these threading technologies, the factor-of-2 reduction will only allow usage of half the physical cores. Again, this provides additional safety that the user's machine will unlikely crash due to workload.

> **Q: What are the recommended system requirements for this GUI?**

A: Most functionality apart from Executor does not require immense resources. Running in a multi-processing manner with the Executor Module, on the other hand, should be used with the following guidelines:
- 2.5 GB of RAM per active core allocated towards a study. **DO NOT** allocate too many cores towards a study/studies if your RAM capacity is insufficient. The program will become slow and unresponsive, as your computer will run out of memory.
- An acceptable cooler for your CPU, as it will be going at it for ~15 minutes per single subject visit at "High" quality setting for both Structural and ASL modules. For a study of 120 visits split among all cores in a 6-core machine, that amounts to ~5 hours of heavy workload. Assuming the user wishes to utilize their machine to the fullest extent without overheating, a recommended top-end air-cooler is the [NH-D15 Chromax Black](https://noctua.at/en/nh-d15-chromax-black) for 6-16 physical core CPUs.
- As file writing/reading occurs in the course of the analysis, having an SSD or NVMe-SSD over a traditional harddrive can allow for reduced processing time.

---

## License

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

- **[GPL license](https://www.gnu.org/licenses/gpl-3.0-standalone.html)**
