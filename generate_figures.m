function [data] = generate_figures()
% generate_figure_main: main script for generating figures for manuscript
% Author: William Davis Haselden
%         The Pennsylvania State University, University Park, PA

close all; clear all; clc;
currentFolder = pwd;
addpath(genpath(currentFolder));
fileparts = strsplit(currentFolder, filesep);
if ismac
    rootfolder = fullfile(filesep, fileparts{1:end},'NOFeedbackData');
else
    rootfolder = fullfile(fileparts{1:end},'NOFeedbackData');
end
addpath(genpath(rootfolder))

%% perform data analysis
[data] = generate_structures(rootfolder);

%% Figure 1
% illustration of Krogh cylinder and input parameters of model (no data)
%% Figure 2
fig_2_a_b(data); %perivascular NO (uniform, regional, proximal) & GC activation curve
fig_2_c_d_e(data); %Vessel diameter, NO production, GC activation plots (uniform, regional, proxiomal)
%% Figure 3
fig_3_a_b_c(data); %perivascular NO, O2, CcO activity (uniform, regional, proximal)
fig_3_d_e_f(data); %Vessel diameter, NO production, CcO inhibition plots (uniform, regional, proxiomal)
%% Figure 4
fig_4_a_b_c_d(); %example time dependent setup: NO production, dilation kernel, GC activation curve, NO sensitvity (m)
fig_4_e_f(data); %example time dependent setup: NO in SM, vasodilation, relative undershoot
%% Figure 5
fig_5_a_b(data); %model dynamics similar to 0.7s and 10 s whisker puff - Drew PNAS 2011 - impulse response undershoot
%% Figure 6

fig_6_a(data); %vessel oscillation continue beyond 6 second kernel when m>>1
fig_6_b(data); %frequency response of sinusoidal NO production
fig_6_c(data); %frequency response of white gaussian noise NO production rate
fig_6_d(data); %Increased NO degradation during vasodilation is a linear system
fig_6_e(data); %Baseline GC activation could alter the power of vasomotion (at 0.16 Hz)

%% Figure 7

fig_7_a(data); %altering NO sensitivity; slope, m = 1/3/5
fig_7_b(data); %altering plasma free Hemoglobin, pfHgb = 1,20,40 uM
fig_7_c(data); %altering Hematocrit, Hct = 20%, 45%, 60%
fig_7_s(data); %static vs dynamic

fig_7_static(data);

end

