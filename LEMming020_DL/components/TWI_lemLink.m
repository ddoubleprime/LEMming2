%% TWI_lemLink.m - LEMming component that generates a timeseries of precipitation regions
%% and returns a grid of rainwater discharge in m/yr the
%% same size as 'topo', which is used in weighting the flow accumulation
%% grid.

%% V.016, for LEMming 016, DW 1 June 2012.

% Remember that "rain rates" here are treated as runoff rates. So if your
% infiltration fraction is 0.5, divde your rain rates by 2. Eventually,
% hopefully, the model will handle infiltration more elegantly

% Here we should either read in the file from a pre-generated timeseries,
% or synthesize the correct spatial gradient.

if READ_IN_PRECIP

    % PrecipGrid would update from a file here, but we're going to use a static pattern
    % so we should load that in the pre-execution code for efficiency
    
else    % synthesize the precipitation history
    
    % Let's start simple, with a step change based on elevation:

    % Assign the accumulated rain rates to the appropriate region
    PrecipGrid(:,:) = Prate1;    % reset the overall precip grid

    % Update inversion level here
    if READ_IN_INVLEVEL
    % For now assumes that the variable 'ilhist' contains the inversion level
    % history relative to initial inversion level in steps identical to the external 
    % integration timestep.
        ilind = round( (tmax-t) / dt_ext ) + 1; % assumes SL history is longer than tmax
        ilind(ilind < 1) = 1; % No zero index; 
        z_inv = z_inv_final + ilhist(ilind); % history should be referenced to final IL
    else
        % Slave to sealevel
        z_inv = sealevel+inv_ht;   % meters in topographic reference frame; this will be controlled by a boundary history
    end
    
    belowInv = topo < z_inv;  % mask of sub-inversion cells
    PrecipGrid(belowInv) = Prate2; 

end