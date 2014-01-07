clear all; close all; clc;
addpath(genpath('../Libraries/eeglab12_0_2_5b'));
f = ls('sets/*.set');

eeglab; clc;

%this for loop iterates over all datasets in dataset1. i've commented out
%the iteration because we will implement the algorithms for one dataset
%first.

for fid = 1%:size(f,1) 
   
   fn = cat(2,strtrim(f(fid,:)));
   EEG = pop_loadset('filename',fn,'filepath','C:\\Users\\Administrator\\Documents\\MATLAB\\bcicomp\\sets\\');
   zhang; %perform zhang's method
end


%%
% classification;