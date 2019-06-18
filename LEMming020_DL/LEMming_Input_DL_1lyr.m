% LEMming - a MATLAB-based landscape evolution model
% Copyright (C) 2012-2015 Dylan Ward
% 
% 
% Developer can be contacted at dylan.ward@uc.edu
% 
% and 
% 
% University of Cincinnati, 
% Dept. of Geology, 
% Cincinnati, OH 45221
% 
% 
% This program is free software; you can redistribute it and/or 
% modify it under the terms of the GNU General Public License as 
% published by the Free Software Foundation; either version 2 of 
% the License, or (at your option) any later version. See License.txt 
% and gpl-2.0.txt for the full terms of this license.
% -------------------------------------------------------------------

% LEMming_Input.m - Input script to set parameters for a LEMming (V020 and up) run.
% DETACHMENT-LIMITED version only (LEMming_DL) - does not track regolith or
% sediment, with the exception of rockfall debris.

clear all
close all
clc

Setup_Name = 'DL';  % As in, LEMming_Input_<Setup_Name>.m
run_name = 'prelim';

% Load a topo file
inFile = 'landscapes/GenericRidge.mat';

% grid size input

x = 240;     % pixels x
y = 200;     % pixels y
dx = 10;      % meters per pixel
dy = dx;


% grid z parameters
z_scale_input = .001;     % Vertically exaggerate the topography by this much.


% Initial max elevation (m). It is useful to specify this here so it can be
% used to define stratigraphy and plotting limits before the input
% topography is loaded. It does not directly affect the topography itself.
z_init = 500;              

% boundary conditions
% Boundary types available: 1) fixed-elevation 2) periodic (BUGGY) 3) follow minimum. If periodic is
% specified for one side of a boundary pair, the setting for the other side
% is ignored

BOUNDARY_NORTH = 1;
BOUNDARY_SOUTH = 1;
BOUNDARY_EAST = 3;
BOUNDARY_WEST = 1;

z_bound = 1;          % Fixed-elevation boundary elevation (m)
borderwidth = 5;      % Width of edge over which boundary condition is applied, cells

OCEAN_BOUND = 0;      % Use variable sea-level boundary condition
READ_IN_SEALEVEL = 0; % Reads a sealevel history from a file (static for now)

% File containing sea level history in variable 'CompositeSL' in intervals
% that are exactly the same as dt_ext and the first index as modern SL
% (zero)
SL_file = 'DeltaSLInvLev';  

% timesteps and model duration
VARIABLE_DT_TOGGLE = 1;     % use a variable timestep
dtDefault = 1e-5;        % initial timestep, yrs; dynamically adjusts
dtMax = 1000;             % max possible timestep, years
dtMin = 1e-2;                  % Minimum timestep (years) - used by all dynamically-stepped processes

dt_safety_factor = 5;   % divides limiting timestep. set lower for faster, less accurate, possibly unstable run (should be stable if > ~4).

z_max_plot = 0;                 % Maximum elevation used for plotting. Set negative or zero to instead use the max elevation of the initial topo.
cinterval = 5;                  % contour interval for plotting, m
tmax = 10e6;                     % years
plottime = 10000;               % years between plots
tracktime = .1 * plottime;      % years between recording the current erosion rates

% shaded relief plot controls
VE = 4;                        % (times) Vertical exaggeration of plot
lightpos = [10 5 10];
ambient = 0.5;
az = -50;
el = 50;

% Toggle saving mode. 
% 0: Only saves beginning and end states.
% 1: Saves each figure, and beginning and end states.
% 2: Saves entire workspace each time a figure is plotted.
SAVEMODE = 2;          

% physical parameters
g = 9.82;                  % m/s^2, gravitational acceleration
f = 0.4;                   % Darcy-Weisbach friction factor. Here we use a value chosen to give ~1 m/s flow for a 50-cm-deep channel at a slope of 0.01. 
rho_w = 1000;              % Density of water, kg/m^3
rho_rock = 2700;           % Density of rock, kg/m^3
rho_reg = 1500;            % Density of regolith/rockfall debris, kg/m^3
    
% Parameters for the channel geometry:
stream_WDR = 10; % dimensionless stream width-to-depth ratio; assumed constant after Finnegan et al. (2005) and Wobus et al. (2006)

% base level control        
rock_uplift = 1e-4;      % Initial/default rock uplift rate relative to boundary condition (m/yr)
VAR_ROCK_UPLIFT = 1;       % toggle to use variable rock uplift rate matrix defined below
% rock uplift hist: t_chg(yr) rate(m/yr). Uses default until first defined
% time point.

% rock_uplift_mtx = [ 1.8e6,  0;
%                     2.0e6,  4.375e-5 ];
                
rock_uplift_mtx = [ 1.8e6,  0;
                    2.0e6,  1e-4;
                    2.1e6,  0;
                    2.5e6,  1e-4;
                    2.6e6,  0;
                    3.0e6,  1e-4;
                    3.1e6,  0 ;
                    3.5e6,  1e-4;
                    3.6e6,  0;
                    3.8e6,  1e-4;
                    3.9e6,  0;
                    4.2e6,  1e-4;
                    4.3e6,  0;
                    4.5e6,  1e-4;
                    4.6e6,  0;
                    4.9e6,  1e-4;
                    5.0e6,  0;
                    5.5e6,  1e-4;
                    5.6e6,  0;
                    6.0e6,  1e-4;
                    6.2e6,  0;
                    7.5e6,  1e-4;
                    7.6e6,  0;
                    8.8e6,  1e-4;
                    9e6,  0   ];
% rock_uplift_mtx = [ 0.2e6,  0;
%                     0.4e6,  1e-4;
%                     0.6e6,  0;
%                     0.8e6,  1e-4;
%                     1.0e6,  0;
%                     1.2e6,  1e-4;
%                     1.4e6,  0;
%                     1.6e6,  1e-4;
%                     1.8e6,  0;
%                     2.0e6,  1e-4;
%                     2.2e6,  0;
%                     2.4e6,  1e-4;
%                     2.6e6,  0;
%                     2.8e6,  1e-4;
%                     3.0e6,  0;
%                     3.2e6,  1e-4;
%                     3.4e6,  0;
%                     3.6e6,  1e-4;
%                     3.8e6,  0;
%                     4.0e6,  1e-4;
%                     4.2e6,  0;
%                     4.4e6,  1e-4;
%                     4.6e6,  0;
%                     4.8e6,  1e-4;
%                     5e6,  0   ];
                
% Sealevel tie points: depending how sea level history is treated, one of
% these values will be used to reference the history to the model domain.
% This should be simplified in the future.
sealevel_final = 1169.6;  % meters elevation of starting SL, topo reference frame
        %% -1169.6 is the original minimum of the extracted/extrapolated
        %% topo/bathymteric swath used to make the idealized Hawaii landscape
sealevel_init = sealevel_final;% - (tmax)*.0026;  % meters elevation of starting SL, topo reference frame


% climate parameters
UNIFORM_PRECIP = 1; % Rainrate is uniform and constant across the whole grid
RRateMean = 1;  % Mean rainrate (m/yr) used only if UNIFORM_PRECIP is true

% TWI_lemLink parameters (step change model), top to bottom rates:
Prate1 = 0.6;   % m/yr runoff (50% loss, 1.2 m/yr mean from 20km - 40km)
Prate2 = 2.7;   % m/yr runoff (50% loss, 5.3 m/yr mean from 0km - 20km)
inv_ht = 1400;  % m relative to sealevel
z_inv_init = inv_ht + sealevel_init;  % initial inversion height, meters elevation, model reference frame. Modern 20km precip break = 1400m
z_inv_final = inv_ht + sealevel_final;

READ_IN_INVLEVEL = 0;
% File containing inversion level history in variable 'IL_var' in intervals
% that are exactly the same as dt_ext and the first index as modern IL
% (zero)
IL_file = 'DeltaSLInvLev'; 
IL_var = 'InvLevel'; % 'CompositeIL'

READ_IN_PRECIP = 0; % Reads a precipitation grid from a file (static for now)
PrecipFile = 'PrecipGrid521';  % Filename containing "PrecipGrid" variable, same size as topo
F_Runoff = 0.5;         % Fraction of precip that results in runoff

%%%%% Rainmaker_lemLink parameters
% % % Precip regions are defined using the same paradigm as stratigraphy. This
% % % list could be generated by a random-events generator according to
% % % distributions in space, time, and intensity, then loaded here.
% %               %= [Pmean   xp1   yp1    xp2   yp2  tpstart tpend]
% % Precip_regions = [2.5     1     y/2+1  x     y      0     tmax  ;...
% %                   0.5     1     1      x     y/2    0     tmax     ];

BEDROCK_FLUV_MODE = 4;  
% rocktype 0 - default substrate
k0 = 1e-3;              % fluvial erodibility constant
m0 = 0.5;               % Area exponent for basic stream power
n0 = 1;                 % Shear stress exponent for detachment-limited erosion
kappa0 = 1e-2;           % Hillslope transport coefficient (m/year) - scale rate. This is the depth-averaged velocity of regolith at a slope of 1 (45 degrees) if transport is linear in slope.
sc0 = tand(30);         % Critical slope (L/L) at which transport rate is infinite (Roering)
rfslope0 = inf;         % Positive slope (L/L) above which qualifies a rockfall source



% additional rocktypes, beginning with 1
%        [index k       m   n kappa sc rfslope]
RTproplist = [1 k0/200 0.5  1 1e-4 tand(60)  1.7; ...      % caprock
              2 k0/200 0.5  1 1e-4 tand(30)  inf; ...      % hard debris
              3 k0     0.5  1 1e-2 tand(30)  inf;          % soft debris
              4 k0/10  0.5  1 5e-4 tand(30)  inf]';        % medium debris

            % Toggle erosion rate summing
            % 1: Sum hillslope and channel components for each cell, proportionally to the fraction of the cell occupied by each
            % 2: Use maximum of hillsope or channel erosion, corrected for fraction of cell that is each
            % 3: Use maximum uncorrected for stream fraction of cell; e.g.,
            %    if stream erosion is higher, entire cell is treated as a channel %%DEFAULT%%
            EROSION_SUM_MODE = 1;   
          
%%% Avalanching - may apply to rockfall debris
AVALANCH_TOGGLE = 1;
angleOfRepose = 30 ;             % degrees
avalanchFreq = .5 ;              % average number per year



% Rockfall parameters
DO_ROCKFALL = 1;

RF_Debris_Rtype = 2;    % Which of the above defined rocktypes represents rockfall debris?
RF_Debris_PlotCutoff = 0.0; % m - thicknesses of debris less than this are plotted as the underlying substrate

% Rockfall source
RFSource_curv = 0; % Curvature below which qualifies a rockfall source (universal - threshold slope is set by rock type.)

DepoAngleCutoff = 25;   % Degrees - angle of repose for rockfall debris

Backwearing_Rate = inf;    % m/yr - Backwearing erosion rate cap on rockfall
RF_Ht = 10; % m - Typical height of a rockfall event, usually the thickness of the hardcap is a good default

% Convert assigned backwearing rate to a volumetric erosion rate
RFerodeRate = Backwearing_Rate * (x * dx) * RF_Ht;

% additional stratigraphy
% [rt x1 y1 z1 x2 y2 z2]
%bookcliffs160 = [1 1 1 1700 160 160 1600; 1 1 1 1100 160 160 900; 1 1 1 700 160 160 650; 1 1 1 550 160 160 400; 1 1 1 300 160 160 200];
%bookcliffssimple100 = [1 1 1 1200 100 100 1100; 1 1 1 0 100 100 -100];
%multicliffs = [1 1 1 1700 x y 1600; 1 1 1 1100 x y 900; 1 1 1 700 x y 650; 1 1 1 550 x y 400; 1 1 1 300 x y 200];
%twocliffs = [1 1 1 -300+75 x y -300+65; 1 1 1 -300+30 x y -300+25];
%smallcliff = [1 1 1 z_init+10 x y z_init-50];
%nocliff = [1 1 1 z_init+1000 x y z_init+999];
%clams = ones(1000,7);clams(:,2) = ceil(rand(1000,1).*x-1);clams(:,3) = ceil(rand(1000,1).*y-1);clams(:,4) = ceil(rand(1000,1).*z_init);clams(:,4) = clams(:,4);clams(:,5) = clams(:,2)+1;clams(:,6) = clams(:,3)+1;clams(:,7) = clams(:,4)-1;clams(clams == 0) = 1;
%tiltstrat = ones(y,7);tiltstrat(:,3) = 1:1:y;tiltstrat(:,4) = .015 * dy .* (tiltstrat(:,3)-round(4*y/5)) + z_init;tiltstrat(:,6) = tiltstrat(:,3);tiltstrat(:,5) = tiltstrat(:,5).*x;tiltstrat(:,7) = tiltstrat(:,4)-50;tiltstrat(tiltstrat == 0) = 1;
%lefttiltstrat = ones(x,7);lefttiltstrat(:,2) = 1:1:x;lefttiltstrat(:,4) = -.03 * dx .* (lefttiltstrat(:,2)-round(3*x)) - 60;lefttiltstrat(:,5) = lefttiltstrat(:,2);lefttiltstrat(:,6) = lefttiltstrat(:,6).*y;lefttiltstrat(:,7) = lefttiltstrat(:,4)-10;lefttiltstrat(lefttiltstrat == 0) = 1;
backtiltstrat = ones(y,7);backtiltstrat(:,1) = 1;backtiltstrat(:,3) = 1:1:y;backtiltstrat(:,4) = -.06 * dy .* (backtiltstrat(:,3)-round(y)) + -300;backtiltstrat(:,6) = backtiltstrat(:,3);backtiltstrat(:,5) = backtiltstrat(:,5).*x;backtiltstrat(:,7) = backtiltstrat(:,4)-20;backtiltstrat(backtiltstrat == 0) = 1;
backtiltstrat2lyr1 = backtiltstrat;
backtiltstrat2lyr2 = backtiltstrat2lyr1;
backtiltstrat2lyr2(:,4) = backtiltstrat2lyr1(:,7) - 30;
backtiltstrat2lyr2(:,7) = backtiltstrat2lyr1(:,7) - 40;

% notch the center by setting the middle column rocktype to a weaker type
%notch = [0 x/2 1 z_init+1000 x/2+1 y z_init-2000];
notch1 = ones(y,7);notch1(:,1) = 4;notch1(:,2) = x/2;notch1(:,3) = 1:1:y;notch1(:,4) = -.06 * dy .* (notch1(:,3)-round(y)) + -300;notch1(:,5) = x/2+1;notch1(:,6) = notch1(:,3);notch1(:,7) = notch1(:,4)-20;
notch2 = notch1;
notch2(:,4) = notch1(:,7) - 30;
notch2(:,7) = notch1(:,7) - 40;

backtiltstrat1lyr = [backtiltstrat; notch1];
backtiltstrat2lyr = [backtiltstrat2lyr1; backtiltstrat2lyr2; notch1; notch2];
stratlist = backtiltstrat1lyr';

%%%%%%%%%%%%%%%%%%%%%%% CALL EXECUTION COMPONENTS AND RUN THE MODEL

% Set path so components can be found
addpath('components','components/upslope','tools','landscapes');

%% CALL ROUTINE TO SET UP THE RUN VARIABLES %%
SetVars_lemLink

%% Now, run the core program %%
LEMming020


%%%%%% End of LEMming_Input.m %%%%%%----------------------------------------