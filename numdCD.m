function out = numdCD(all_samplexrf) 
% Function that calculates the number density for elements of interest for
% each sample. Outputs to an array
%
% Written by Allie Koester and Nat Lifton, Purdue University, Feb 2022
% koestea@purdue.edu
%Copyright 2022, Purdue University
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License, version 2,
% as published by the Free Software Foundation (www.fsf.org).

% Input must be major element oxides in weight percent for each sample
% in the following sequence
% SiO2, TiO2, Al2O3, Fe2O3, FeO, MnO, MgO, CaO, Na2O, K2O, and P2O5

%check to make sure input data is correct before proceeding 
n=11;
if length(all_samplexrf) < n
    error('Error \nNumber of columns must be %d',n); %display error message
elseif length(all_samplexrf) > n
    error('Error \nNumber of columns must be %d',n);
end 

%check if calcite/carbonate was submitted - assume leftover weight percent is CO2
s = sum(all_samplexrf,2); %sum up each sample weight percent
num_samples = height(s);
m = 98; %arbitrary minimum total weight percent for silicates (should add to ca. 100%) 
for i = 1:num_samples
    if s(i) < 60 %maximum summed weight percent for carbonates
        all_samplexrf(i,12) = 100-s(i); 
    elseif s(i) < m
        error('Error \nTotal summed weight percents for each sample must be over %d',m); %display error message and exit
    else
        all_samplexrf(i,12) = 0; %add a placeholder for non-carbonate samples
    end
end

all_samplexrf = all_samplexrf./100; %convert to fractional value
% num_samples = height(all_samplexrf);
ND = zeros(num_samples,12); 

all_sample_SiO2 = all_samplexrf(:,1); 
all_sample_TiO2 = all_samplexrf(:,2); 
all_sample_Al2O3 = all_samplexrf(:,3);
all_sample_Fe2O3 = all_samplexrf(:,4);
all_sample_FeO = all_samplexrf(:,5);
all_sample_MnO = all_samplexrf(:,6);
all_sample_MgO = all_samplexrf(:,7);
all_sample_CaO = all_samplexrf(:,8);
all_sample_NaO = all_samplexrf(:,9);
all_sample_K2O = all_samplexrf(:,10);
all_sample_P2O5 = all_samplexrf(:,11);
all_sample_CO2 = all_samplexrf(:,12); 

%Use these elemental fractions for measured oxide values from XRF data
load Elfrac.mat
%first column of Elfrac structure is fraction of oxygen, second column is fraction of element
%in the oxides (top to bottom) of SiO2, TiO2, Al2O3, Fe2O3, FeO, MnO, MgO, CaO, Na2O,
%K2O, P2O5, CO2


NA = 6.02214076E23; %Avogadro's number

%Sum up the weight percents, calculate normalization to 1.0 for each sample
for a =1:num_samples
    norm = 1./sum(all_samplexrf,2); %input file, start with data 
end 


%calculate number density for each element
%need to separate out from oxygen

ND(:,1) = (Elfrac(1,2).*all_sample_SiO2.*NA)./28.085; %Si number density (ND)
ND(:,2) = (Elfrac(2,2).*all_sample_TiO2.*NA)./47.8670; %Ti ND
ND(:,3) = (Elfrac(3,2).*all_sample_Al2O3.*NA)./26.9815; %Al ND 
ND(:,4) = (Elfrac(4,2).*all_sample_Fe2O3.*NA)./55.845; %Fe2O3 ND
ND(:,5) = (Elfrac(5,2).*all_sample_FeO.*NA)./55.845; %FeO ND
ND(:,6) = (Elfrac(6,2).*all_sample_MnO.*NA)./54.938004; %Mn ND
ND(:,7) = (Elfrac(7,2).*all_sample_MgO.*NA)./24.3055; %Mg ND
ND(:,8) = (Elfrac(8,2).*all_sample_CaO.*NA)./40.078; %Ca ND
ND(:,9) = (Elfrac(9,2).*all_sample_NaO.*NA)./22.9898; %Na ND
ND(:,10) = (Elfrac(10,2).*all_sample_K2O.*NA)./39.0983; %K ND
ND(:,11) = (Elfrac(11,2).*all_sample_P2O5.*NA)./30.973762; %P ND
%Add up all oxygen fractions
ND(:,12) = (Elfrac(1,1).*all_sample_SiO2 + Elfrac(2,1).*all_sample_TiO2 +...
    Elfrac(3,1).*all_sample_Al2O3 + Elfrac(4,1).*all_sample_Fe2O3 +...
    Elfrac(5,1).*all_sample_FeO + Elfrac(6,1).*all_sample_MnO +...
    Elfrac(7,1).*all_sample_MgO + Elfrac(8,1).*all_sample_CaO +...
    Elfrac(9,1).*all_sample_NaO + Elfrac(10,1).*all_sample_K2O +...
    Elfrac(11,1).*all_sample_P2O5 + Elfrac(12,1).*all_sample_CO2).*NA./15.9994;

%Normalize results to 100% and save in an array
%format: Si Ti Al Fe1 Fe2 Mn Mg Ca Na K P O
results = ND.*norm; 

out = results;
end 

