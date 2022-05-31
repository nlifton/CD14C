function out = LSDscalingCD(h,Rc,SPhi,w,consts,ND) 

% Implements the Lifton Sato et al scaling scheme for spallation.
%
% Syntax: scalingfactor = LSDscalingCD(h,Rc,SPhi,w,consts,ND);
%
% Where:
%   h = atmospheric pressure (hPa)
%   Rc = cutoff rigidity (GV)
%   SPhi = solar modulation potntial (Phi, see source paper)
%   w = fractional water content of ground (nondimensional)
%   ND = number density of each element per measured major element oxide
%   composition
%   
% Vectorized. Send in scalars or vectors of common length. 
%
% Written by Nat Lifton 2011, Purdue University
% Based on code by Greg Balco -- UW Cosmogenic Nuclide Lab
% balcs@u.washington.edu
% April, 2007
% Part of the CRONUS-Earth online calculators: 
%      http://hess.ess.washington.edu/math

% Copyright 2001-2007, University of Washington
% All rights reserved
% Developed in part with funding from the National Science Foundation.

% Modified by Allie Koester and Nat Lifton, 2022 Purdue University
% koestea@purdue.edu
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License, version 2,
% as published by the Free Software Foundation (www.fsf.org).

% convert pressure to atmospheric depth

X = h.*1.019716;

CRef = consts.P14nRef_q + consts.P14pRef_q; %qtz theoretical spallation PR 
% at SLHL reference, i.e., h = 1013.25, Rc = 0, SPhi = 624.5718,w = 0.066.

EthRef = consts.ethfluxRef;
ThRef = consts.thfluxRef;

% Full version with cross sections and includes thermal and
% epithermal fluxes

NSite = NeutronsCD(h,Rc,SPhi,w,consts,ND); 
PSite = ProtonsCD(h,Rc,SPhi,consts,ND);
[ethflux,thflux] = NeutronsLowE(h,Rc,SPhi,w);

%Nuclide-specific scaling factors as f(Rc) corresponding to time vector in
%ScalingLSD_CD

% Compositionally dependent scaling factor for sample composition relative to 100% SiO2
% composition, referenced to SLHL
Site.Ccd = (NSite.P14CD + PSite.P14CD)./CRef; %part of Equation 3


Site.E = NSite.E;%Nucleon flux energy bins
Site.eth = ethflux./EthRef; %Epithermal neutron flux scaling factor as f(Rc)
Site.th = thflux./ThRef;%Thermal neutron flux scaling factor as f(Rc)

out = Site;

