function out = make_consts_CD14C()

% This function creates and saves a structure with relevant constants and
% external data for the C and Be exposure age and erosion rate calculators.  
%
% Syntax: make_conts_CD14C
% (no arguments)
%
% Written by Greg Balco -- Berkeley
% Geochronology Center
% balcs@u.washington.edu -- balcs@bgc.org
% 
% Modified by Brent Goehring -- Lamont-Doherty Earth Observatory
% goehring@ldeo.columbia.edu
% and Nat Lifton -- Purdue University
% nlifton@purdue.edu
% April, 2011
% %
% Copyright 2011, University of Washington, Columbia University and the
% Purdue University
% All rights reserved
% Developed in part with funding from the National Science Foundation.

% Modified by Allie Koester and Nat Lifton -- Purdue University 2022
% koestera@purdue.edu

% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License, version 2,
% as published by the Free Software Foundation (www.fsf.org).
%

consts.version = '3.4_CD14C - 05/2022'; 
consts.prepdate = fix(clock);

% C-14 decay constant 

consts.l14 = log(2)/5700;

% Note that the uncertainty is not used in exposure-age or erosion-rate
% calculators. Here only for development purposes. 

consts.dell14 = log(2) * 40/(5700^2);

% Effective attenuation length for spallation in rock
% Commonly accepted value: 160 g/cm2
% For discussion see Gosse and Phillips (2001)

consts.Lsp = 160;

% Reference production rates at SLHL for spallation according to LSDn
% scaling framework (Lifton et al., 2014)
%

% CRONUS Primary Results from Borchers et al. (2016) and Phillips et al.
% (2016) and Young et al. (2014)
% Note that correct SPhi values per Eq. 3 of Lifton et al. (2014) are not included in Borchers
% et al. (2016) or Marrero et al. (2016) CRONUScalc code, but global means of CRONUS primary 
% site production rates (calculated individially by Lifton, 2/2016) differ only in 2nd 
% decimal place from calculation using same method (individual site PRs, then global mean) 
% but with incorrect SPhi values. 
% LS uses trajectory-traced Rc values for <14ka from Pavon-Carrasco et al. (2014) in Lifton (2016); 
% LD uses geocentric dipole from Pavon-Carrasco in Lifton (2016). 
% Sample concentrations recalculated using Hippe and Lifton (2014) - 
% Straight sample means calculated first, then site production rates.
% Values are straight means of all site PRs and standard deviations - NL 2/22

consts.P14_ref_LS = 13.50;  %nuclide-specific LSD - Trajectory Traced Rc - PC14
consts.delP14_ref_LS = 0.89; %nuclide-specific LSD - Trajectory Traced Rc - PC14
consts.P14_ref_LD = 13.71;  %nuclide-specific LSD - Geocentric dipole - PC14
consts.delP14_ref_LD = 1.20; %nuclide-specific LSD - Geocentric dipole - PC14

% Atomic number densities (atoms target/g mineral)
% (Moles target * Avogadro's Number / Formula weight)

% Quartz (SiO2)

consts.NatomsQtzO = 2.00458052253686E+22; %adjusted by AJK/NL - uses 2019 exact definition of Avogadro's number
consts.NatomsQtzSi = 1.00229026126843E+22; %adjusted by AJK/NL - uses 2019 exact definition of Avogadro's number

% Spallogenic Nuclide Production Cross-Sections (n & p)
% References:
%   Reedy 2013; Reedy 2007, LPSC 38, abstract 1192; Reedy, 2011, personal communication
%   Jull et al., 1998
%   Imamura et al., 1990
%   JENDL - Fukahori et al., 2002; Watanabe et al., 2011
%   TENDL - Koning et al., 2019

load react14C % compiled by AJK
% Use variable names as indicated
consts.E = react14C.E;
consts.O16nn2pC14 = react14C.O16nn2pC14R; %Reedy 2013; Imamura et al., 1990
consts.O16pxC14 = react14C.O16pxC14R; %Reedy 2013; Reedy 2007, LPSC 38, abstract 1192; Reedy, 2011, personal communication
consts.Si28nxC14 = react14C.Si28nxC14R; %Reedy 2013; Reedy 2007, LPSC 38, abstract 1192; Reedy, 2011, personal communication
consts.Si28pxC14 = react14C.Si28pxC14R; %Reedy 2013; Reedy 2007, LPSC 38, abstract 1192; Reedy, 2011, personal communication
consts.Al27nxC14 = react14C.Al27nxC14R; %Jull et al., 1998; Reedy 2007, LPSC 38, abstract 1192; Reedy, 2011, personal communication
consts.Al27pxC14 = react14C.Al27pxC14R; %Jull et al., 1998; Reedy 2007, LPSC 38, abstract 1192; Reedy, 2011, personal communication
consts.Ca40nxC14 = react14C.Ca40nxC14J; %JENDL
consts.Ca40pxC14 = react14C.Ca40pxC14J; %JENDL
consts.Fe56nxC14 = react14C.Fe56nxC14R; %Jull et al., 1998; Reedy 2007, LPSC 38, abstract 1192; Reedy, 2011, personal communication
consts.Fe56pxC14 = react14C.Fe56pxC14R; %Jull et al., 1998; Reedy 2007, LPSC 38, abstract 1192; Reedy, 2011, personal communication
consts.K39nxC14 = react14C.K39nxC14J;%JENDL
consts.K39pxC14 = react14C.K39pxC14J;%JENDL
consts.Mg24nxC14 = react14C.Mg24nxC14R; %%Jull et al., 1998; Reedy 2007, LPSC 38, abstract 1192; Reedy, 2011, personal communication
consts.Mg24pxC14 = react14C.Mg24pxC14R; %%Jull et al., 1998; Reedy 2007, LPSC 38, abstract 1192; Reedy, 2011, personal communication
consts.Mn55nxC14 = react14C.Mn55nxC14J; %JENDL
consts.Mn55pxC14 = react14C.Mn55pxC14J;%JENDL
consts.Na23nxC14 = react14C.Na23nxC14J;%JENDL
consts.Na23pxC14 = react14C.Na23pxC14J;%JENDL
consts.P31nxC14 = react14C.P31nxC14J;%JENDL
consts.P31pxC14 = react14C.P31pxC14J;%JENDL
consts.Ti48nxC14 = react14C.Ti48nxC14J;%JENDL
consts.Ti48pxC14 = react14C.Ti48pxC14J;%JENDL
consts.P31nxC14T = react14C.P31nxC14T; %TENDL-2019 
consts.P31pxC14T = react14C.P31pxC14T; %TENDL-2019 
consts.Na23nxC14TJ = react14C.Na23nxC14TJ; %spliced data Low E TENDL to 200 MeV, High E JENDL - to compare JENDL & TENDL
consts.Na23pxC14TJ = react14C.Na23pxC14TJ; %spliced data Low E TENDL to 200 MeV, High E JENDL - to compare JENDL & TENDL

% Paleomagnetic records for use in time-dependent production rate schemes
% Derived from Nat Lifton's compilation of paleomagnetic data from
% various sources. See Lifton et al. (2014) and Lifton (2016)

load PavonPmag 
%Pavon-Carrasco et al. (2014) to 14kyr, GLOPIS-75 to 75 ka, PADM2M to 2Ma, 1950-2010 Rc  
% From Lifton, 2016

% Relative dipole moment and time vector - Pavon-Carrasco et al. (2014)
consts.GLOPADMM0 = GLOPADMM0; 
consts.tM = tM; 

consts.PavonRc = PavonRc; % data block 
consts.lat_Rc = lat_Rc; % lat and lon indices for Rc data block
consts.lon_Rc = lon_Rc;
consts.tRc = tRc; % time vector for Rc data block

% 
% Effective pole positions and field strengths inferred from Pavon-Carrasco et al. (2014) field
% reconstructions for last 14000 yr. They are for
% the same times as the RC slices in the data block above. 
% Lifton, 2016

consts.MM0_Pavon = MM0_Pavon;
consts.lat_pp_Pavon = lat_pp_Pavon;
consts.lon_pp_Pavon = lon_pp_Pavon;

% Solar variability

% Solar variability from Usoskin et al. 2011
% 0-11400 yr - 100-yr spacing, 0 yr is 2010

% Per Lifton et al. (2014)
%Convert Usoskin et al. (2011) solar modulation potential to Sato Force Field Potential due to different
%assumed Local Interstellar Spectrum and other factors

SPhi = 1.1381076.*SPhi - 1.2738468e-4.*SPhi.^2;

consts.SPhi = SPhi;
consts.SPhiInf = mean(SPhi);

%Calculated reference PRs using hi-resolution cross-sections (800 energy
%bins instead of original 200, per Argento et al., 2015) 7/20/15
load ReferencePFluxCD %AJK updated P14nRef_q and P14pRef_q values 
%Reference values for scaling via Sato et al. (2008) spectra

consts.P14nRef_q = P14nRef_q; %14C neutron reference production in SiO2
consts.P14pRef_q = P14pRef_q; %14C proton reference production in SiO2

consts.nfluxRef = nfluxRef; %integral reference neutron flux
consts.pfluxRef = pfluxRef; %integral reference proton flux
consts.ethfluxRef = ethfluxRef; %integral reference epithermal flux
consts.thfluxRef = thfluxRef; %integral reference thermal neutron flux
consts.mfluxRef = mfluxRef; % reference muon flux components

% Finish up

save("consts_CD14C",'consts');

disp(['Constants version ' consts.version]);
disp('Saved'); 


