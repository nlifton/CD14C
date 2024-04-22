function out = CD14C(sampledata)
% wrapper script to calculate the compositionally dependent in situ 14C production rate for a sample
%
% syntax = CD14C('sampledata.txt')
% output is a binary .mat file

% Written by Allie Koester and Nat Lifton 2022, Purdue University
% koestea@purdue.edu

% Based on code by Greg Balco (Balco et al., 2008)
% April, 2007
% Part of the CRONUS-Earth online calculators: 
%      http://hess.ess.washington.edu/math
% Copyright 2001-2007, University of Washington
% All rights reserved
% Developed in part with funding from the National Science Foundation.

% Also based on code by Nat Lifton (Lifton et al., 2014; Lifton, 2016);
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License, version 2,
% as published by the Free Software Foundation (www.fsf.org).

version = '1.0 - 05/2022';

%% Import data

%Input data format is Sample_Name; Latitude (deg N); Longitude (deg E); Elevation (m); 
% SiO2, TiO2, Al2O3, Fe2O3, FeO, MnO, MgO, CaO, Na2O, K2O, and P2O5 all in
% weight percent

FID = fopen(sampledata);
data = textscan(FID,'%s %n %n %n %n %n %n %n %n %n %n %n %n %n %n');
fclose(FID);
dstring='';

%Import the sample data vectors
all_sample_name = data{1}; %sample name
all_lat = data{2}; %latitude
all_long = data{3}; %longitude
all_elv = data{4}; %elevation
all_xrf = cell2mat(data(5:end)); %matrix of all the sample major element data

%calculate elemental number density for each sample 
all_ND = numdCD(all_xrf);

num_samples = length(all_lat);

%% Load the constants file
% make_consts_CD14C;
load consts_CD14C.mat;

%% Create a sample structure for each sample 
for a = 1:num_samples
    sample.sample_name = all_sample_name{a};
    sample.lat = all_lat(a);
    sample.long = all_long(a);
    sample.elv = all_elv(a);
    sample.ND = all_ND(a,:);

    sample.pressure = ERA40atm(sample.lat,sample.long,sample.elv);

%% Scaling and Production calculations 

    disp(sample.sample_name); %display current sample name
    scaling14 = ScalingLSD_CD(sample,consts);
    
    %Calculate the spallation production rate using gridded Rc values
    P14_CD(a,:) = scaling14.SF_LS_CD.*consts.P14_ref_LS;  % Equation 3
                                    %LSD calibrated qtz spallation PR

    %Calculate the spallation production rate using geocentric dipole
    %(similar to UW v3 online calculator approach (Balco et al., 2008))
    P14gd_CD(a,:) = scaling14.SF_LD_CD.*consts.P14_ref_LD; % Equation 3
                                      %LSD calibrated qtz spallation PR
     

end

%% output the results
%save output data to workspace
out.CDVersion = version;
out.ConstsVersion = consts.version;
out.ScalingVersion = scaling14.version;
out.ID = all_sample_name';
out.tv = scaling14.tv; %time vector, relative to 2010 (t = 0)
out.P14_CD = P14_CD; %time- and compositionally dependent site production 
% rate relative to that geologically calibrated for quartz
out.P14gd_CD = P14gd_CD; %time- and compositionally dependent site production 
% rate relative to that geologically calibrated for quartz, for geocentric dipole


save(replace(sampledata,".txt","results.mat"),'out');

disp(['Constants version ' consts.version]);
disp('Saved'); 

end