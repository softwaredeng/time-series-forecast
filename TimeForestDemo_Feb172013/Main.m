% code is NOT optimized
% the results may be different from each run due to randomness
% run mexC.m first (only need to run once)


clc;clear;close all;
[allErr,allForest]=TForest();
allErr
TemporalImportanceCurve(allForest);
