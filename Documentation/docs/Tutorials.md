
# Tutorials & FAQ

----
## Beginner Tutorials

----
### How to set up an ASL-BIDS dataset


----
### How to run ExploreASL using Matlab

The first thing you have to do, to use **ExploreASL**, is to clone the **ExploreASL** repository. If you want to run **ExploreASL** from Matlab, we recommend to clone the main repository directly from the official [GitHub website](https://github.com/ExploreASL/ExploreASL). You also have the option to download the zipped version or to download an older [release](https://github.com/ExploreASL/ExploreASL/releases).

If you are new to Matlab, we recommend checking out a [Matlab tutorial](https://www.mathworks.com/support/learn-with-matlab-tutorials.html). It can be helpful to add the **ExploreASL** directory to your Matlab paths. Open Matlab, select the **Home** tab, and add the **ExploreASL** directory including its subfolders using the **Set Path** option. Now change your working directory, using the **Current Folder** tab or the **cd** command, to the **ExploreASL** directory.

To run **ExploreASL** you have to type in the following command in the **Command Window**: `ExploreASL`. If you already created an **ASL-BIDS dataset**, you can run the full default **ExploreASL** pipeline like this:

```
DataParPath = 'C:\...\MY-BIDS-DATASET\DataParFile.json';
ProcessData = true;
SkipPause = true;
[x] = ExploreASL_Master(DataParPath, ProcessData, SkipPause);
```

----
### How to run ExploreASL using the GUI


----
### How to run ExploreASL using the docker image



----
## Advanced Tutorials





----
## FAQ



