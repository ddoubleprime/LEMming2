% Enforce_Borders_lemLink.m - applies the spatial boundary conditions and
% tracks border modifications to the mass flux

%% DETACHMENT-LIMITED v020


switch BOUNDARY_NORTH
    case 1

        % Store unmodified N border
        tempN = topo(end-(borderwidth-1):end,:);

        % Flatten N border
        topo(end-(borderwidth-1):end,:) = z_bound; %min(topo(end - borderwidth,:)); %% NORTH
        % Record volumes removed by flattening the borders. Should equal the volume flux into the borders.
        dVN = dVN + sum(sum(tempN)) * CellArea;  % North volume removed by BC

        
    case 2
        
        % BOUNDARY_SOUTH = 2; Handle wraparound in S-code
        
    case 3
        
        % Follow N border - TEMPED TO MAX 7/9/15
        topo(end,:) = max(topo(end-2*borderwidth:end-1,:),[],1); % NORTH
        
        % A better idea - constant slope, tiny gradient down to border
        %topo(end,:) = (topo(end-1,:) + topo(end-2,:)) / 2; % NORTH
        

        
end

switch BOUNDARY_SOUTH
    case 1

        % Store unmodified S border
        tempS = topo(1:borderwidth,:);

        % Flatten S border
        topo(1:borderwidth,:) = z_bound;  % SOUTH
        % Record volumes removed by flattening the borders. Should equal the volume flux into the borders.
        dVS = dVS + sum(sum(tempS)) * CellArea;  % South volume removed by BC


                
    case 2
        
        % Wraparound edges
        % Edge row/col corresponds to second row/col on opposite side

        % NORTH-SOUTH
        topo(1:2,:) = topo(end-1:end,:);

        %StreamWs(1:2,:) = StreamWs(end-1:end,:);
        slopes(1:2,:) = slopes(end-1:end,:);

        AccGrid(end,:) = AccGrid(2,:);
        AccGrid(1,:) = AccGrid(end-1,:);
        
        R(end,:) = R(2,:);
        R(1,:) = R(end-1,:);

        
    case 3
        
        % Follow S border
        topo(1,:) = min(topo(2:2*borderwidth,:),[],1);  % SOUTH

        % gradient
        %topo(1,:) = (topo(2,:) + topo(3,:)) / 2;  % SOUTH
        
        

        
end
        
        

switch BOUNDARY_EAST
    case 1

        % Store unmodified E border
        tempE = topo(:,end-(borderwidth-1):end);
        
        % Flatten E border
        topo(:,end-(borderwidth-1):end) = z_bound; % EAST
        % Record volumes removed by flattening the borders. Should equal the volume flux into the borders.
        dVE = dVE + sum(sum(tempE)) * CellArea;  % East volume removed by BC




    case 2
        
        %BOUNDARY_WEST = 2;     % Handle wraparound in W-code
        
    case 3
 
        % Follow E border
        topo(:,end) = min(topo(:,end - 2*borderwidth:end-1),[],2); % EAST
        %topo(:,end) = mean(topo(:,borderwidth+1:end - borderwidth),2); % EAST
        
        % gradient
        %topo(:,end) = (topo(:,end - 1) + topo(:,end - 2)) / 2; % EAST
        

        
end

switch BOUNDARY_WEST
    case 1

        % Store unmodified W border
        tempW = topo(:,1:borderwidth);

        % Flatten W border
        topo(:,1:borderwidth) = z_bound;  % WEST
        % Record volumes removed by flattening the borders. Should equal the volume flux into the borders.
        dVW = dVW + sum(sum(tempW)) * CellArea;  % West volume removed by BC   


        
    case 2
        
        % Wraparound edges
        % Edge row/col corresponds to second row/col on opposite side

        % EAST-WEST
        topo(:,1:2) = topo(:,end-1:end);

        %StreamWs(:,1:2) = StreamWs(:,end-1:end);
        slopes(:,1:2) = slopes(:,end-1:end);

        AccGrid(:,end) = AccGrid(:,2);
        AccGrid(:,1) = AccGrid(:,end-1); 
        
        R(:,end) = R(:,2);
        R(:,1) = R(:,end-1); 
        
        
    case 3
        
        % Follow W border
        topo(:,1) = min(topo(:,2:2*borderwidth),[],2);  % WEST
        %topo(:,1) = mean(topo(:,borderwidth+1:end-borderwidth),2);  % WEST
        
        % gradient
        %topo(:,1) = (topo(:,2) + topo(:,3)) / 2;  % WEST

        % Apply regolith boundary condition 

        
end
        
        % Don't let anything erode past base level.
        topo(topo < z_bound) = z_bound;     

        % Don't use the border area in calculating the timestep
        dtMask = BorderGrid;
        if OCEAN_BOUND
            dtMask = dtMask + isOcean;
        end
        
