%% Rainmaker.m - standalone development version of random storm generator
%% for LEMming. Based on the stochastic rockfall-cosmogenic nuclide
%% simulator written by DW in 2008. The lemLink version will return a grid
%% of precipitation in m/yr the same size as topo, which is used in
%% weighting the flow accumulation grid.

%% V.016, for LEMming 016, DW 28 Mar. 2011.


% Start doing stuff

if ~UNIFORM_PRECIP

% We'll put a time-list manager here to turn zones on and off as time goes
% on. For now let them all run the whole time.

PrecipGrid(:,:) = 0;    % reset the overall precip grid

for rf_zone = 1:size(Precip_regions,1)

    
    % Read in the properties matrix. Let's get space working first then
    % handle time
    Pmean = Precip_regions(rf_zone,1);    % Desired overall precipitation rate in meters/yr
    Alpha = Precip_regions(rf_zone,2);    % Storm area PDF scaling exponent
    A_min = Precip_regions(rf_zone,3);    % Minimum storm area
    A_max = Precip_regions(rf_zone,4);    % Maximum storm area
    
    xp1 = min(Precip_regions(rf_zone,5),Precip_regions(rf_zone,7)); % ll corner
    xp2 = max(Precip_regions(rf_zone,5),Precip_regions(rf_zone,7)); % ur corner
    yp1 = min(Precip_regions(rf_zone,6),Precip_regions(rf_zone,8)); % ll corner
    yp2 = max(Precip_regions(rf_zone,6),Precip_regions(rf_zone,8)); % ur corner

    % trapping to correct for input coordinates that are off the grid
    xp1(xp1 < 1) = 1;
    xp2(xp2 > x) = x;
    yp1(yp1 < 1) = 1;
    yp2(yp2 > y) = y;
    
    px = (1 + (xp2-xp1)); % x dimension of region
    py = (1 + (yp2-yp1)); % y dimension of region
    Asim = (dx*dy) * (px*py);    % Square meters area of interest
    rr_field = zeros(py,px);     % Make the subgrid
    
% Generate events given a particular scaling law and erosion rate

%% this needs to be rederived to yield a quasi-steady mean precipitation
%% rate over time, independent of timestep. Should look up appropriate
%% distributions to use.
if Alpha == 1.5; Alpha = 1.4999; end
Kappa = ((3-2*Alpha)) / (2*Alpha*(A_min^Alpha)*(A_max^(3/2 - Alpha) - A_min^(3/2 - Alpha)));


% Do events
n_events = round(dt*Asim*Kappa);
    
Pas = rand(1,n_events);
As = A_min.*(1-Pas).^(-1/Alpha);
rr = rand(size(As)) * (Pmean);                   % This is the rainfall rate over the area of this event. Will vary by a PDF. As will duration.


% Now, rain this timestep's events out on the landscape.

Ws = sqrt(As);  % Width of each event when expressed as a square
W_cellsX = round(Ws ./ dx);   % Widths in pixels
W_cellsY = min(W_cellsX,py);
W_cellsX = min(W_cellsX,px);    % Trap for events bigger than the grid

reXs = ceil( (rand(1,n_events) .* (px-1)) + 0.5);    % Event X origin
reYs = ceil( (rand(1,n_events) .* (py-1)) + 0.5);    % Event Y origin

for event = 1:n_events
    
    event_rows = reYs(event):reYs(event)+W_cellsY(event);
        event_rows(event_rows > py) = event_rows(event_rows > py) - py; % Wraparound boundary
    event_cols = reXs(event):reXs(event)+W_cellsX(event);
        event_cols(event_cols > px) = event_cols(event_cols > px) - px; % Wraparound boundary
    
    rr_field(event_rows,event_cols) = rr_field(event_rows,event_cols) + rr(event);
end

% Assign the accumulated rain rates to the appropriate region
PrecipGrid(yp1:yp2,xp1:xp2) = PrecipGrid(yp1:yp2,xp1:xp2) + rr_field; 

end % for rf_zone

end % if ~UNIFORM_PRECIP
