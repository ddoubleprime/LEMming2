%% Rainmaker_lemLink.m - Takes the timeseries of recipitation regions from
%% the input file and returns a grid of rainwater discharge in m/yr the
%% same size as topo, which is used in weighting the flow accumulation
%% grid.

%% V.016, for LEMming 016, DW 13 May 2011.


% Start doing stuff

if ~UNIFORM_PRECIP

PrecipGrid(:,:) = 0;    % reset the overall precip grid

active_zones = find(Precip_regions(:,6) <= t & Precip_regions(:,7) >= t); % Get row indices of the active zones

if ~isempty(active_zones)   % trap for no active rain events

for rfzi = 1:length(active_zones)
    rf_zone = active_zones(rfzi);
    
    % Read in the properties matrix. Let's get space working first then
    % handle time
    Pmean = Precip_regions(rf_zone,1);    % Desired overall precipitation rate in meters/yr
    
    xp1 = min(Precip_regions(rf_zone,2),Precip_regions(rf_zone,4)); % ll corner
    xp2 = max(Precip_regions(rf_zone,2),Precip_regions(rf_zone,4)); % ur corner
    yp1 = min(Precip_regions(rf_zone,3),Precip_regions(rf_zone,5)); % ll corner
    yp2 = max(Precip_regions(rf_zone,3),Precip_regions(rf_zone,5)); % ur corner

    % trapping to correct for input coordinates that are off the grid
    xp1(xp1 < 1) = 1;
    xp2(xp2 > x) = x;
    yp1(yp1 < 1) = 1;
    yp2(yp2 > y) = y;

    
% Assign the accumulated rain rates to the appropriate region
PrecipGrid(yp1:yp2,xp1:xp2) = PrecipGrid(yp1:yp2,xp1:xp2) + Pmean; 

end % for rf_zone

end % if ~isempty(rf_zone)

end % if ~UNIFORM_PRECIP