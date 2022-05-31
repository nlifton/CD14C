function scaling = ScalingLSD_CD(sites,consts)
% 
% scaling is a structure with many fields:
%%
% results.main_version: version of this function
% results.tv: time vector against which to plot Rc and P
%
% LS = Lifton et al. 2014 nuclide-dependent scaling framework, with Bob Reedy Cross-Section Data (9/2010)
%
% results.P_XX: P(t) at site in scaling scheme XX (atoms/g/yr) (vector)

%
% Written by Greg Balco -- UW Cosmogenic Nuclide Lab and Berkeley
% Geochronology Center
% balcs@u.washington.edu
% 
% Modified by Brent Goehring -- Lamont-Doherty Earth Observatory
%    goehring@ldeo.columbia.edu
% and Nat Lifton -- Purdue University
%    nlifton@purdue.edu
% September, 2010

% Modified by Allie Koester and Nat Lifton
%
% Copyright 2022, University of Washington, Columbia University and
% Purdue University
% All rights reserved
% Developed in part with funding from the National Science Foundation.


% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License, version 2,
% as published by the Free Software Foundation (www.fsf.org).

%
scaling.version = '1.0'; 
sfa = ['LS';'LD'];

num_samples = length(sites.lat);
% Make the time vector

% Age Relative to t0=2010
tv = [0:10:50 60:100:3060 3260:200:75060 76060:1000:800060 802060:2000:2000060 logspace(log10(2002060),7,200)];

LSRc = zeros(num_samples,length(tv));
LDRc = zeros(num_samples,length(tv));

% Need solar modulation for Lifton SF's

this_SPhi = zeros(size(tv)) + consts.SPhiInf; % Solar modulation potential for Sato et al. (2008)
this_SPhi(1:length(consts.SPhi)) = consts.SPhi; % Solar modulation potential for Sato et al. (2008)
w = 0.066; % water content for Sato & Niita (2006), gravimetric, about 14% volumetric per Fred Phillips

% interpolate an M for tv > 0 for dipole field. 10/1/15 Pavon et al. SHA.DIF.14k
temp_M1 = interp1(consts.tRc,consts.MM0_Pavon,tv(1:91));
temp_M2 = interp1(consts.tM,consts.GLOPADMM0,tv(92:end));
temp_M = cat(2,temp_M1,temp_M2);

for a = 1:num_samples
    clear sample;
    % Load the input data structure
%     sample.site = sites.site(a);
%     sample.sample_name = sites.site(a);
    sample.lat = sites.lat(a);
    sample.long = sites.long(a);
    sample.elv = sites.elv(a);
    sample.ND = sites.ND(a:end); %number densities for each sample
%    sample.age = sites.truet(a);
    
    % Pressure correction
    
    sample.pressure = ERA40atm(sample.lat,sample.long,sample.elv);

    % Catch confusion with pressure submission. If sample.pressure is already 
    % set, it should have a submitted value. If zero, something is wrong. 
    % This should never happen in online use. 

    if sample.pressure == 0
        error(['Sample.pressure = 0 on sample ' sample.sample_name]);
    end
    
    % catch for negative longitudes before Rc interpolation
    if sample.long < 0; sample.long = sample.long + 360;end


% TIME-VARYING GAD FIELD
%   Lifton et al. 2014 -  Fit to Trajectory-traced GAD dipole field as f(M/M0), as long-term average.
    GDRc = [-448.004 1189.18 -1152.15 522.061 -103.241 6.89901 0];

    % For PavonRc - 10/1/15

    [longi,lati,tvi] = meshgrid(sample.long,sample.lat,tv(1:91)); %Interpolate gridded Rc values per Lifton (2016)
    LSRc(a,1:91) = interp3(consts.lon_Rc,consts.lat_Rc,consts.tRc,consts.PavonRc,longi,lati,tvi);
    LSRc(a,92:end) = temp_M(92:end).*polyval(GDRc,cos(d2r(sample.lat)));

   zIndex = LSRc(a,:) < 0;
   LSRc(a,zIndex) = 0;

   LDLat = abs(90-angdist(sample.lat,sample.long,consts.lat_pp_Pavon,consts.lon_pp_Pavon));
   LDRc(a,1:91) = temp_M(1:91).*polyval(GDRc,cos(d2r(LDLat)));
   LDRc(a,92:end) = temp_M(92:end).*polyval(GDRc,cos(d2r(sample.lat)));

   zIndex = LDRc(a,:) < 0;
   LDRc(a,zIndex) = 0;

    LSD = LSDscalingCD(sample.pressure,LSRc(a,:),this_SPhi,w,consts,sample.ND); %gridded Rc values
    LSDgd = LSDscalingCD(sample.pressure,LDRc(a,:),this_SPhi,w,consts,sample.ND); %geocentric dipole

    SF_LS_CD(a,:) = LSD.Ccd; %scaling factor for gridded Rc values
    SF_LD_CD(a,:) = LSDgd.Ccd; %scaling factor for geocentric dipole


% 5. Results structure assignment
        
end

% Time template
scaling.tv = tv;
scaling.LSD = LSD; %output the elemental production with energies
scaling.LSDgd = LSDgd; %output the elemental production with energies - geocentric dipole
scaling.SF_LS_CD = SF_LS_CD;
scaling.SF_LD_CD = SF_LD_CD;

