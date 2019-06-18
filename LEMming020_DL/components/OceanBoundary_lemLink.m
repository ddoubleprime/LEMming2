%% OceanBoundary_lemLink
%% Enforces a variable sea-level boundary condition, wherein topography
%% below sea level is preserved when sea level rises but all fluxes of
%% sediment are discarded when reaching the sea.

% This stuff is set up once at the beginning (relocated to setVars_lemLink)

% sealevel = sealevel_init;   % meters in topographic reference frame; this will be controlled by a boundary history
% isOcean = topo < sealevel;  % mask of ocean cells
% topo_master = topo;         % intialize topo buffer

% The below gets run every loop after erosion and such

topo(isOcean) = topo_master(isOcean) - Regolith_H(isOcean);     % update working topography below sea level
Regolith_H(isOcean) = 0;    % Strip all regolith entering the ocean

% Update sea level here
if READ_IN_SEALEVEL
    % For now assumes that the variable 'CompositeSL' contains the sea level
    % history relative to initial sea level in steps identical to the external 
    % integration timestep.
    slind = round( (tmax-t) / dt_ext ) + 1; % assumes SL history is longer than tmax
    slind(slind < 1) = 1; % No zero index; 
    sealevel = sealevel_final + CompositeSL(slind); % history should be referenced to final SL
    
else    % synthesize the sea level history
    sealevel = sealevel_init;% + t*.0026;% + 60 * sin(t / 10000);   % meters in topographic reference frame; this will be controlled by a boundary history
end

isOcean = topo < sealevel;  % mask of ocean cells
topo_master(~isOcean) = topo(~isOcean);   % update master topography above sea level
topo(isOcean) = sealevel;   % set new base level