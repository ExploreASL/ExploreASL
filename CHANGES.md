# ExploreASL V1.1.3

----

## Bug Fixes

* hotfix minor bug in running the import using DCMTK without the Matlab Image Processing Toolbox #30

# ExploreASL V1.1.2

----

## Bug Fixes

* hotfix minor bug in loading NIfTIs containing lesion masks in CAT12 #28

----

# ExploreASL V1.1.1

----

## Bug Fixes

* hotfix minor bug in creating participants.tsv #23

----


# ExploreASL V1.1.0

  

----

  

Additional information about ExploreASL can be found on the [ExploreASL](www.ExploreASL.org) website.
To start ExploreASL from Matlab, type in: ```ExploreASL_Master```. 
The walkthrough document and how-to videos can be found on the [](https://sites.google.com/view/exploreasl) [ExploreASL](www.ExploreASL.org) website, and in the [Neuroimage paper]([https://pubmed.ncbi.nlm.nih.gov/32526385/](https://pubmed.ncbi.nlm.nih.gov/32526385/)). Further documentation is work in progress. For any help please contact the main authors/developers at h.j.mutsaerts@amsterdamumc.nl or j.petr@hzdr.de.

  

----

## Bug Fixes

* Bug fixes and overall code improvements related to the BIDS import workflow (#11)
* Registration with poor CBF contrast will not issue an error anymore but correctly switch to control-T1w registration only (#17)
* Unexisting x.Sequence field fixed, now an appropriate warning is issued and this field is defined automatically by ```xASL_adm_DefineASLSequence.m``` (#16)
----

## Features
* Quantification can now be fully disabled by: ```x.ApplyQuantification = [0 0 0 0 0];``` (#14)
* Insert option to disable M0-ASL registration (#13)
* Update of CAT (Computational Anatomy Toolbox) from version 12.5 to 12.7 (#2)

---
## Work in progress
* Minor improvements of custom scripts for BBB-ASL and BIDS (#8)
* Minor improvements regarding unit testing of ExploreASL (#10)
* Additional warnings for ExploreASL users (#12)
----

## Documentation
 
* Recent changes include the improvement of the documentation within the ExploreASL structure using markdown files and the introduction of a new documentation repository (#7)
* Some function headers were added for increased understandability (#19). These can be viewed in Matlab by: ```help ExploreASL_Master``` where you can replace ExploreASL_Master by the actual function name


----

# ExploreASL V1.0.0

  ----
This is the first release version.
