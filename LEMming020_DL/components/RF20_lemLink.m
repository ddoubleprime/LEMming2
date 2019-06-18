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

% Stochastic rockfall runout simulator. This version gets run inline with
% LEMming v020 and above.

RFdz(:,:) = 0;
RFH_exp(:,:) = 0;
RFH = RF_Debris_H;      % temp RF array for the AOR code

curves = del2(topo);

sourceMask = curves < RFSource_curv & slopes > rfslope;

[sourceYs,sourceXs] = find(sourceMask);
sourceZs = topo(sourceMask);

n_sources = length(sourceZs);
n_events_do = n_sources;

eventList = ceil( rand(1,n_events_do) .* (n_sources-.1) );
[eventListFilt, uniqueEventIdx] = unique(eventList);  % Screen for dupes

erodeBucket = erodeMaxRate * dt;    % Figure out max erosion allowed in this timestep



% Now simulate the collection of events
if ~isempty(eventList);
                    %temprary:
               % disp('Rockfall!')
        
    for event = uniqueEventIdx'
        if erodeVol < erodeBucket;
            
            
        sourceX = sourceXs(event);
        sourceY = sourceYs(event);
        sourceZ = sourceZs(event);

        % Erosion. Maintain relief; subtract source cell to height of highest cell
        % below it.
        event_rows = sourceY-1:sourceY+1;
            event_rows(event_rows > y) = event_rows(event_rows > y) - y; % Wraparound boundary
            event_rows(event_rows < 1) = event_rows(event_rows < 1) + y; % Wraparound boundary
        event_cols = sourceX-1:sourceX+1;
            event_cols(event_cols > x) = event_cols(event_cols > x) - x; % Wraparound boundary
            event_cols(event_cols < 1) = event_cols(event_cols < 1) + x; % Wraparound boundary
        
        localDrops = sourceZ - topo(event_rows,event_cols);
        erodeHt = min(min( localDrops(localDrops > 0) )); % could use mean of these.
        if isempty(erodeHt); erodeHt = 0; end
        

        % Rockfall should expand here to its lower density. However, this
        % is somewhat problematic - if the topography expands with it, it
        % will cause rockfall debris to spread over the lip of the cliff;
        % if not, then the thicknesses are not handled consistently with
        % the gradient-based angle of repose code. Will need to solve, perhaps with
        % a flag for new rockfall - or an expansion based on the delivery
        % of zero-age material once the age-tracking code is working.
         RFH(sourceY,sourceX) = erodeHt; %* (rho_rock/rho_reg);
         RFH_exp(sourceY,sourceX) = erodeHt * (rho_rock/rho_reg - 1);
         RFAge(sourceY,sourceX) = t;    % record model time of event

        
        % Volume tracking for regulation of erosion rate
        eventVol = erodeHt * CellArea;
        erodeVol = erodeVol + eventVol;
        
        else break % the event loop
        end  % if erodebucket
        
    end  % for event
        
    
    % Deposition - use the shallow landsliding routine to distribute
    % debris at angle of repose
        LSMODE = 1;         % zero for regolith, one for rockfall  - sets mode of the shallow LS code
        ShallowLS_lemLink
        
        % Update erosion/deposition amount grid: layer handler likes to have the difference. 
        RFdz = RFH - RF_Debris_H + RFH_exp;  % material due to expansion is added here and distributed next timestep.

end % if ~isempty(eventList)

% No accumulation in the border areas.
RFdz(logical(BorderGrid)) = 0;

erodeVol = erodeVol - erodeBucket; % get excess erosion
erodeVol = max(erodeVol,0); % reset to zero if no leftovers