
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
### How to run a compiled ExploreASL Version


To compile ExploreASL you have to run the `xASL_adm_MakeStandalone.m` script. If necessary, you can also ask the developer team for a specific compiled version. Providing a compiled version for every operating system and corresponding Matlab version is currently not feasible for us. Please feel free to ask us for help though!

#### Windows

Let's assume you want to run the compiled version of **ExploreASL v1.4.0**. Check the contents of the folder created by `xASL_adm_MakeStandalone.m`, which contains the compiled version. There should be a file called `xASL_1_4_0.exe`. We recommend using the command line interface now. For this you can go to the address bar of your file explorer. Type in `cmd` to open the command prompt in the current folder. The following command will start **ExploreASL** and process the dataset of the DataParFile:

```
xASL_1_4_0.exe "c:\path_to_your\DataParFile.json"
```

The executable will extract all necessary data from the CTF archive within the folder. This is totally normal.

#### Linux

On Linux you can basically do the same as above. Let's assume we chose the option to create a compiled **ExploreASL** version labelled with the **latest** tag within `xASL_adm_MakeStandalone.m` script. We can run the ExploreASL shell script with a specified Matlab MCR (here we use version 96 e.g.) using the following command:

```
/bin/bash /path/to/xasl/xASL_latest/xASL_latest.sh /path/to/mcr/v96/ "path_to_your_DataParFile.json"
```


----
### How to run ExploreASL using the docker image

First you have to pull an official docker image from the **ExploreASL** repository:

```
docker pull docker.pkg.github.com/exploreasl/docker/xasl:1.x.x
```

If you want to rename the docker image, you can use the docker tag command:

```
docker tag docker.pkg.github.com/exploreasl/docker/xasl:1.x.x xasl:1.x.x
```

To start a docker container of a dataset you can use the following command:

```
docker run -e DATAPARFILE="DataSet/DataParFile.json" -v /home/.../incoming:/data/incoming -v /home/.../outgoing:/data/outgoing xasl:1.x.x
```

- Here DATAPARFILE is an environment variable which is a relative path to the data parameter file of your dataset.
- ```/home/.../incoming:/data/incoming``` is used to mount your dataset folder (```/home/.../incoming```) to its corresponding docker dataset folder (```/data/incoming```). 
- The same notation is used to mount the docker dataset output folder (```/data/outgoing```) to its corresponding real output folder on your drive (```/home/.../outgoing```).


----
## Advanced Tutorials





----
## FAQ



