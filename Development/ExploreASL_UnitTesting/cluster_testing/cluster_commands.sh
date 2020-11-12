# basic test showing matlab version
# module load matlab/R2020a
# matlab -batch "v=version;display(v);exit;"


# download xASL 
git clone https://github.com/ExploreASL/ExploreASL.git

# copy relevant data
cp ./ExploreASL/External/TestDataSet TestDataSet

# submit job
sbatch xASL_job.run
