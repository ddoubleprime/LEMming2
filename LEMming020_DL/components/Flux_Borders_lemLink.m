% Flux_Borders_lemLink.m - applies the flux boundary conditions and
% tracks border modifications to the mass flux

defbndflx = 0;%max(max((Qout)))+1; % Just enough to be divergent at base level

switch BOUNDARY_NORTH
    case 1

        Qout(end,:) = defbndflx;

    case 2
        
        BOUNDARY_SOUTH = 2;     % Handle wraparound in S-code
        
end

switch BOUNDARY_SOUTH
    case 1

        Qout(1,:) = defbndflx;
                
    case 2
        
        % Periodic N-S
        Qout(end,:) = Qout(2,:);
        Qout(1,:) = Qout(end-1,:);
        
end
        
        

switch BOUNDARY_EAST
    case 1

        Qout(:,end) = defbndflx; 

    case 2
        
        BOUNDARY_WEST = 2;     % Handle wraparound in W-code
        
end

switch BOUNDARY_WEST
    case 1

        Qout(:,1) = defbndflx; 
        
    case 2

        % Periodic E-W
        Qout(:,end) = Qout(:,2);
        Qout(:,1) = Qout(:,end-1); 

end

