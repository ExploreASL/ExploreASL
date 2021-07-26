clear all;
close all;
clc;
v=version;
disp(['Matlab version is ', v]);
cd ExploreASL/
x=ExploreASL_Master('../TestDataSet/analysis/DataParameters_LowQ.json',true,true);
