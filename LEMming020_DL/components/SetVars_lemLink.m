% LEMming - a MATLAB-based landscape evolution model
% Copyright (C) 2011 Dylan Ward
% 
% 
% Developer can be contacted at djward@unm.edu 
% 
% and 
% 
% University of New Mexico, 
% Dept. of Earth and Planetary Sciences, 
% Albuquerque, NM 87131
% 
% 
% This program is free software; you can redistribute it and/or 
% modify it under the terms of the GNU General Public License as 
% published by the Free Software Foundation; either version 2 of 
% the License, or (at your option) any later version. See License.txt 
% and gpl-2.0.txt for the full terms of this license.
% -------------------------------------------------------------------

%% A MATLAB-based, heavily vectorized LEM with the ability to add layers of
%% variable lithology

%% D. Ward, 12/16/2007: Versions 002-015
%% D. Ward, 5/25/2011: Version 016
%% Version 016 is a major upgrade. Improvements include:
%%      - Significant improvements to regolith routing and transport,
%%      bedload, and suspended sediment fluxes
%%      - Bedrock channel erosion now works as expected
%%      - Spatially and temporally variable rainfall rate
%%      - Now uses a nested timestepping scheme whereby each process
%%      determines its own optimal timestep, and the external integration
%%      timestep is fixed.
%%      Now includes MEX-components written in C as replacements for two
%%      slow loops in the flow-routing and regolith-routing procedures.
%%      - Initial topography must now be generated externally and provided
%%      in a MAT-file; the inline topo generator was becoming cumbersome to
%%      maintain.


% Variables should be specified in LEMming_Input_<setupname>.m
% or otherwise set by script before running LEMming

% DETACHMENT-LIMITED VERSION -LEMming020-DL


%%%%%%%%%%%%%%%%%%%% EXECUTION %%%%%%%%%%%%%%%%%%%%

% Prompt for a model run name
run_name = input('What would you like to call this model run? > ','s');
disp('Creating folder for this model run...')
run_filename = ['RUNS/' Setup_Name '/' run_name ' ' datestr(now,'dd-mmmm-yyyy HH.MM')];
mkdir(run_filename);        % Create a subfolder for this run

% Copy current codes into the new subfolder so the run can be replicated.
% NEEDS MAD DEBUGGING - COPY ALL COMPONENTS FOR INSTANCE
copyfile('LEMming020.m',['./' run_filename]);
% copyfile('tools/MakeLEMmingMov.m',['./' run_filename]);
% copyfile('components/RF_Spread_lemLink.m',['./' run_filename]);
copyfile(['LEMming_Input_' Setup_Name '.m'],['./' run_filename]);

disp('Executing...')

rand('twister', sum(100*clock));  %#ok<RAND> % Reseed the random number generator


% Make XY absolute position grids
    [Xs, Ys] = meshgrid(1:x,1:y);
    Xs = Xs .* dx;
    Ys = Ys .* dy;

% Grid reference areas
    CellArea = (dx * dy);
    GridArea = (x * y) * CellArea;
    GridDelta = mean(dx,dy);  % Use as "universal" grid spacing when dx ~= dy
    
    

if exist('inFile') %#ok<EXIST>
    
    load(inFile);
    %topo = topo';
    topo = topo(1:y,1:x); % Subset to X,Y.
    topo = fill_sinks(topo);    % Fill DEM sinks
    topo = topo * z_scale_input;    % Vertical stretch

   %dx = cellsize; dy = dx; % Use cellsize from file (comment out to override to dx specified above)
    
   
    topo = topo - min(min(topo)) + z_bound;  % Assert that base level is z_bound, as specified above. Usually zero.

else

    disp 'Please specify a valid input landscape MAT-file, and make sure it is on the LEMming path.'
    
end % if inFile




% Dynamic timestep filter
topoPad = zeros(y+2,x+2);
dtfilt = ones(3) / 9 ;

% Avalanching stuff
if AVALANCH_TOGGLE || DO_ROCKFALL
    [rws,cls] = size(topo) ;
    dHReposeReg = dx*tan(angleOfRepose*pi/180) ;
    dHReposeRF = dx*tan(DepoAngleCutoff*pi/180) ;

    dZidx_down = zeros(rws,cls);
    dZidx_up = zeros(rws,cls);
    dZidy_left = zeros(rws,cls);
    dZidy_right = zeros(rws,cls);

    delHdn = zeros(rws,cls) ; delHup = zeros(rws,cls) ;
    delHlt = zeros(rws,cls) ; delHrt = zeros(rws,cls) ;
    
    RFAge   = zeros(rws,cls);
    RFAgeup = zeros(rws,cls);
    RFAgedn = zeros(rws,cls);
    RFAgelt = zeros(rws,cls);
    RFAgert = zeros(rws,cls);
    
end

dZdx = zeros(y,x+1);  % Will create zero-slope margins, could be changed as needed
dZdy = zeros(y+1,x);


% initial conditions

% Miscellaneous constants
    e = exp(1);
    spyr = 60*60*24*365;    % Seconds per year
    
    
% Prepare the stratigraphy grids
    stratlist = [stratlist; 1:length(stratlist(1,:))];   %add unit ids to strat list row 8
    numStypes = length(RTproplist(1,:));

% Initialize property grids
    k = ones(size(topo)) .* k0;          
    m = ones(size(topo)) .* m0;           
    n = ones(size(topo)) .* n0;        
    kappa = ones(size(topo)) .* kappa0;          
    sc = ones(size(topo)) .* sc0;          
    rfslope = ones(size(topo)) .* rfslope0; 

    RTGrid = zeros(size(topo));
    rtype_here = zeros(size(topo));
   

    % Need staggered grids for Kappa and sc. Probably should be
    % populated by the stratigraphy handler
    kappa_x = kappa0 * ones(size(dZdx));
    kappa_y = kappa0 * ones(size(dZdy));
    sc_x = sc0 * ones(size(dZdx));
    sc_y = sc0 * ones(size(dZdy));
    sThresh_x = sc_x-0.00001*sc_x;   % Prevents infinite erosion rates and backdiffusion at critical slope. 
    sThresh_y = sc_y-0.00001*sc_y;   % Prevents infinite erosion rates and backdiffusion at critical slope. 
    
% Initialize the surface stratigraphy layers
    RF_Debris_H = zeros(size(topo));
    RFH = zeros(size(topo));  
    RFH_exp = zeros(size(topo));
    
    dzCum = zeros(size(topo));
    RFdz = zeros(size(topo));

    BRtopo = topo-RF_Debris_H;
    
% Initialize flow-related grids
    if READ_IN_PRECIP && ~UNIFORM_PRECIP
        load(PrecipFile)
        PrecipGrid = PrecipGrid*F_Runoff;   % Scale for infiltration/runoff
    else
        PrecipGrid = ones(size(topo)) * RRateMean;  % Default is uniform rainrate
    end
    
    if READ_IN_INVLEVEL
            load(IL_file, IL_var);
            ilhist = eval(IL_var);
            clear(IL_var);
            ilind = round( tmax / dt_ext ) + 1; % assumes SL history is longer than tmax
            ilind(ilind < 1) = 1; % No zero index; 
            z_inv = z_inv_final + ilhist(ilind); % history should be referenced to final SL
    else
            z_inv = z_inv_init;
    end
    
    slopes = zeros(size(topo));
    AccGrid = zeros(size(topo));
    StreamWs = ones(size(topo));
    BRChan = false(size(topo));
    U = zeros(y,x) + rock_uplift;       % Generic rock uplift rate field
    
% Initialize the time loop and counters
    t = 0;
    dt = dtDefault;
    t_plot = 0;
    t_track = 0;
    erodeVol = 0;
    trackstep = 0;
    NoDistTargets = 0;
    
% Check boundary options and prepare boundaries
    
    % Periodic boundaries have to match the opposite side
    if BOUNDARY_NORTH == 2 || BOUNDARY_SOUTH == 2
        BOUNDARY_NORTH = 2;
        BOUNDARY_SOUTH = 2;
    end
    if BOUNDARY_EAST == 2 || BOUNDARY_WEST == 2
        BOUNDARY_EAST = 2;
        BOUNDARY_WEST = 2;
    end  

    BorderGrid = zeros(size(topo));
    
    if BOUNDARY_NORTH == 1
        BorderGrid(end-(borderwidth-1):end,:) = 1;  % NORTH
    elseif BOUNDARY_NORTH == 3
        BorderGrid(end,:) = 1;
    end
    if BOUNDARY_SOUTH == 1
        BorderGrid(1:borderwidth,:) = 1;    % SOUTH
    elseif BOUNDARY_SOUTH == 3
        BorderGrid(1,:) = 1;
    end
    if BOUNDARY_EAST == 1
        BorderGrid(:,end-(borderwidth-1):end) = 1;  % EAST
    elseif BOUNDARY_EAST == 3
        BorderGrid(:,end) = 1;
    end
    if BOUNDARY_WEST == 1
        BorderGrid(:,1:borderwidth) = 1;    % WEST
    elseif BOUNDARY_WEST == 3
        BorderGrid(:,1) = 1;
    end

    BorderArea = CellArea * sum(sum(BorderGrid));
    BorderGrid = logical(BorderGrid);
    
    % Enforce fixed boundaries
    topo(BorderGrid) = z_bound;
    

    % Variable in time rock uplift rate
    if VAR_ROCK_UPLIFT && ~isempty(rock_uplift_mtx);
        ruindx = 0;
        if t >= rock_uplift_mtx(ruindx+1,1)
           ruindx = ruindx + 1;
           rock_uplift = rock_uplift_mtx(ruindx,2);
           U(:,:) = rock_uplift;       % Generic rock uplift rate field
        end
    else
        VAR_ROCK_UPLIFT = 0;    % disable variable rock uplift if empty matrix
        % uses default rock_uplift rate defined above        
    end
        
    
    % Ocean Boundary Condition
    if OCEAN_BOUND
        if READ_IN_SEALEVEL
            load(SL_file,'CompositeSL')
            slind = round( tmax / dt_ext ) + 1; % assumes SL history is longer than tmax
            slind(slind < 1) = 1; % No zero index; 
            sealevel = sealevel_final + CompositeSL(slind); % history should be referenced to final SL
        else
            sealevel = sealevel_init;
        end
        
        isOcean = topo < sealevel;  % mask of ocean cells
        topo_master = topo;         % intialize topo buffer
        topo(isOcean) = sealevel;   % set new base level
    end

% Set up erosion rate tracking

    dVS = 0; %#ok<NASGU>
    dVN = 0; %#ok<NASGU>
    dVE = 0; %#ok<NASGU>
    dVW = 0; %#ok<NASGU>
    
    Streamgrid_lemLink;             % Establish the flow matrices

    %Enforce_Borders_lemLink;    % Apply boundary modifications to the initial grids
    topo_init = topo;       % For double-checking erosion volumes later
    
    dVS = 0;
    dVN = 0;
    dVE = 0;
    dVW = 0;
    
    numRecSteps = ceil(tmax/tracktime) + 1;
    FLUVtrackVol = zeros(1,numRecSteps);    % Fluvial bedrock erosion volume
    REGtrackVol = zeros(1,numRecSteps);     % Volume of regolith produced
    REGDZtrackVol = zeros(2,numRecSteps);   % Change in elevation due to regolith motion
    RFtrackVol = zeros(2,numRecSteps);      % Change in elevation due to rockfall activity
    TIMEtrackVol = zeros(1,numRecSteps);    % Time of each measurement
    
    
    clear numRecSteps;

    if SAVEMODE == 1 || SAVEMODE == 2
        stateNo = 0;
    end

% Set up slope filter
    window = 3;
    filtpatch = ones(window) / window^2;

% Miscellaneous plot preparation    

    if VAR_ROCK_UPLIFT
        z_max_plot = max(z_max_plot,max(max(topo))) + sum([rock_uplift_mtx(1,1); diff(rock_uplift_mtx(:,1)); tmax-rock_uplift_mtx(end,1)] .* [rock_uplift; rock_uplift_mtx(:,2)]);
    else
        z_max_plot = z_init;
        %z_max_plot = max(z_max_plot,max(max(topo))) + rock_uplift * tmax;
    end
    
    
    %colordef none
    
%% END OF SetVars_lemLink.m %%    