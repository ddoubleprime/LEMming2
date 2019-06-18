% REGOLITH_lemLink.m
% Component for LEMming to handle regolith production and transport by
% hillslope and fluvial processes





            %%%%%%%%% THIS IS THE REGOLITH HANDLER CODE %%%%%%%%%
tStop = t + dt; % This finishes the internal loop and returns to the external program
tReg = t;       % Set the internal time to the current external time   
exposed_time(:,:) = 0;  % reset the grid of bedrock exposure time
dzRegCum(:,:) = 0;   % reset cumulative elevation change

while tReg < tStop;

        %Enforce_Borders_lemLink  % Enforce boundary conditions
        Streamgrid_lemLink       % Update the flow matrices
            
        % Mobile regolith thickness
        Mobile_H = min(Regolith_H,Max_Mobile_H);
        %Mobile_H = Max_Mobile_H;
        %Mobile_H = Regolith_H;
      
        %%% Regolith transport %%%
        
        % Regolith transport rate, Roering 2001
        % First limit slopes to 0.5 degree less than critical so as to
        % prevent infinite transport
        sLimited = slopes;
        sLimited(slopes >= sThresh) = sThresh;
        UregDiff = (kappa .* slopes) ./ ( 1-( (sLimited./sc).^2) );
        %UregDiff = (kappa .* slopes) ./ ( 1-( (slopes./sc).^2) );

            % Note that because of the kappa*slopes term, rate still
            % increases linearly with slope above the critical threshold
            % defined by sThresh

        % Calculate diffusive fluxes
        QregDiff = UregDiff .* Mobile_H .* dx;
        
        % Regolith transport rate by stream transport - old stream power
        % equation
%         UregFluv = kt .* AccGrid.^mt .* abs(slopes).^nt;    

        % Regolith transport rate by stream transport - shear stress
        % equation (Meyer-Peter & Mueller)
        QregFluv = (kt .* StreamWs .* (max(0,StreamTau - TauCrit50)).^nt);    % Volumetric transport rate

        % SuspendedLoad_lemLink     % Calculate suspended and wash fluxes
        % (NOT CURRENTLY FUNCTIONAL)
        QregSusp=0;
        
        % Find bedrock channels where fluvial processes dominate transport
        % and no immobile regolith/sediment is present
        %%% NOTE: is reg_H criterion needed or wanted? 
        %%% NOTE2: Could this be moved to Streamgrid component?
        %%% NOTE3: Area threshold makes sure both fluxes are defined
        %%%     meaningfully, but may imact divide structure?
        QregPotential=UregDiff * Max_Mobile_H * dx;
        %QregPotential=QregDiff;
        if DETACH_LIMIT_STREAMS == 4;
            BRChan = ( (QregFluv + QregSusp > FluxRatioDetach * QregPotential) & AccGrid > 2*dx*dy);% & Regolith_H <= Max_Mobile_H );
        end
        
            
        % Total regolith flux out of cells
        Qout = QregFluv + QregSusp + QregDiff;  % volume flux of sediment out 
        
        
        Flux_Borders_lemLink  % Enforce flux at boundaries
        %% Against my will I'm going to use a loop here and see how slow
        %% things get at the expense of a really cool regolith routing
        %% procedure. This should be stable enough for pretty big timesteps if it works.
        %% Maybe I'll find a way to vectorize it in the future.
        
%         Qin(:,:) = 0;
%         Qadd(:,:) = 0;
% 
%         
%         for yi = 2:y-1
%             for xi = 2:x-1
% 
%               if Qout(yi,xi) > 0    % only do this for cells with an outflux
%                   
%                 QareaY = yi-1:yi+1;  % Only need the immediate neighbors
%                 QareaX = xi-1:xi+1;
% 
%                 pixel_numpatch = pixel_number(QareaY,QareaX);
%                 from = pixel_number(yi,xi);
%                 weightpatch = Ts(pixel_numpatch,from);
%                 weightpatch = -reshape( weightpatch,3,3);
% 
% 
%                 % Multiply the fractional influences by the flux out of
%                 % yi,xi
%                 Qadd = weightpatch * Qout(from);    
%                 Qadd(2,2) = 0;    % Don't add back the outflux from the cell
%                 
%                 Qin(pixel_numpatch) = Qin(pixel_numpatch) + Qadd;   %#ok<SAGROW> % Accumulate influxes
%                 
%               end   % if Qout
%                 
%             end % for xi
%         end % for yi

        % The C version: literally thousands of times faster.
        Qin = LocalQ(x,y,Qout-QregSusp,Ts);

        % Under periodic boundary conditions only, the influx grid needs to
        % be modified before calculating erosion rates
            if BOUNDARY_NORTH == 2 || BOUNDARY_SOUTH == 2
                    % Periodic N-S influx
                        Qin(end,:) = Qin(2,:);
                        Qin(1,:) = Qin(end-1,:);
            end

            if BOUNDARY_EAST == 2 || BOUNDARY_WEST == 2
                    % Periodic E-W influx
                        Qin(:,end) = Qin(:,2);
                        Qin(:,1) = Qin(:,end-1);
            end            

        
        % Sum contributions to the flux divergence
        %Qout(isOcean) = inf;
        Edot_reg = (Qin - Qout) ./ (dx*dy);
        Edot_reg(isOcean) = 0;
        Edot_reg(BorderGrid) = 0;
        Edot_reg(BRChan) = 0;

        %% Figure out regolith production rates for inclusion in timestep
        %% consideration
        ProduceReg_lemLink
        
            

% LEMming uses process-specific timestepping, such that regolith transport, fluvial
% erosion, rockfall, etc. can each have their own optimal timesteps,
% integrated at the external LEMming timestep. These timesteps cannot be
% larger than the external timestep. These timesteps default to the
% external timestep if they cannot be calculated or if VARIABLE_DT_TOGGLE
% is set to zero.

% Here is the internal regolith timestep calculation.

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% CALCULATE TIMESTEP
        
        if ( VARIABLE_DT_TOGGLE )
            
            %% now that we know the rate of change in surface heights due to  
            %% regolith motion and production we need to know over 
            %% what period of time we can project forward with these rates and 
            %% maintain stability of the land surface.  The basic idea here is that
            %% we don't want to take a timestep any longer then it would take to 
            %% reverse the land surface slope between two cells, such that regolith 
            %% should be flowing in the other direction.  In fact, let's make our 
            %% timestep much less than that.
            
            %% this calculation sets the timestep such that the change
            %% in surface elevation nowhere exceeds a set fraction
            %% of the local standard deviation in ice surface elevations
            
            % include regolith production
            dHdtTot = Edot_reg;
            adHdt = abs(dHdtTot) ;
            
            % something like standard deviation of 3x3 cell areas around each cell
            topoPad(2:end-1,2:end-1) = topo;
            topoMean = filter2( dtfilt, topoPad, 'valid' ) ;
            dHmax = sqrt( filter2( dtfilt, (topoMean - topo).^2 ) ) ;
            
            % find limiting timestep for each considered cell
            ind = ( adHdt~=0 & dHmax~=0 & ~BorderGrid) ;    % exclude boundary areas!
            dtLimits = dHmax(ind)./adHdt(ind) ;
            [dtH, idt] = min( dtLimits ) ;
            
            
            % Regolith thickness should not exceed a preset amount in a
            % timestep or the production function calculation becomes
            % inaccurate
            
            dtR = min(min(Reg_prod_limit ./ (rdot_H(rdot_H ~= 0) + rdot_dilate(rdot_H ~= 0))));
            
            % limit timestep to the stripping time or some fraction of the calculated timestep
            if ~isempty(dtR)
                dtReg = min(dtH/2,dtR);
            else
                dtReg = dtH/2;
            end
            
            % Now, we need to make sure that a) the timestep is not larger
            % than the external timestep, and b) that the cumulative internal time
            % doesn't overshoot the external timestep.
            
            if dtReg < dt_min
                
                warndt = sprintf('Warning: REGinfinity desired timestep is smaller\n than the user-specified minimum.\n Instability may result.') %#ok<NOPTS>
                dtReg %#ok<NOPTS>
                dt_min 
                dtReg = dt_min; % Don't go smaller than the specified limit

                
            end
            
            dtReg = min(dtReg,dt);     % Don't go bigger than the external timestep
            
            if (dtReg + tReg) > tStop
                dtReg = tStop - tReg;
            end
            
            % catch an error
            if isempty(dtReg)  
                dtReg = dt ;
            end
            
        else

            dtReg = dt ;
                        
        end %- VARIABLE_DT_TOGGLE -%
            
        % Update the internal time
        tReg = tReg + dtReg;
            
%%% Here's where we solve regolith flux iteratively, so that more regolith
%%% is not transferred out of a cell than the sum of what is in the cell
%%% and the current influx. Because the outfluxes affect the influxes, we
%%% must iterate to converge on a solution within a tolerance specified by
%%% the user.
% go_flag = true;
% loopCount = 0;
% activeCount = activeCount+1;
% while go_flag
%     
%     loopCount = loopCount+1;
%     Qin0 = Qin;   % Store previous influx value for comparison    

% Subdivide timestep so regolith erosion rates are not applied to
% bedrock   
            nzEdot = Edot_reg ~= 0; % nonzero erosion rates mask
            dt_prime = zeros(y,x);  % time it takes to strip regolith
            dt_prime(nzEdot) = (Regolith_H(nzEdot)) ./ -Edot_reg(nzEdot);  % time needed to strip regolith
            dt_prime(dt_prime > dtReg) = dtReg;       % if not stripped, use normal dt
            dt_prime(dt_prime < 0) = dtReg;         % if negative, then thickness increases, use normal dt
            dt_prime(~nzEdot & Regolith_H > 0) = dtReg;        % if no reg erosion, no bedrock erosion
            
        % Update outflux grid
%         Qout = Qout .* (dt_prime./dtReg) + Qin .* ( (dtReg-dt_prime)./dtReg );
%         
%         %Flux_Borders_lemLink  % Enforce flux at boundaries
%         
%         % Recalculate the influx grid using the updated outfluxes
%         Qin = LocalQ(x,y,Qout-QregSusp .* (dt_prime./dtReg),Ts);
% 
%         % Under periodic boundary conditions only, the influx grid needs to
%         % be modified before calculating erosion rates
%             if BOUNDARY_NORTH == 2 || BOUNDARY_SOUTH == 2
%                     % Periodic N-S influx
%                         Qin(end,:) = Qin(2,:);
%                         Qin(1,:) = Qin(end-1,:);
% 
%             end
% 
%             if BOUNDARY_EAST == 2 || BOUNDARY_WEST == 2
%                     % Periodic E-W influx
%                         Qin(:,end) = Qin(:,2);
%                         Qin(:,1) = Qin(:,end-1);
% 
%             end            
% 
%         
%         % Sum contributions to the flux divergence
%         Edot_reg = (Qin - Qout) ./ (dx*dy);
%         Edot_reg(isOcean) = 0;  
%         Edot_reg(BorderGrid) = 0;
%         Edot_reg(BRChan) = 0;
%         
%         QNZ = Qin0 ~= 0;
%         % Test to see if we've converged
%         MaxdQ = max(max(abs( (Qin0(QNZ)-Qin(QNZ)) ./ Qin0(QNZ) )));
%         if isempty(MaxdQ) || MaxdQ < Flux_Tolerance || loopCount >= loopLimit;
%             go_flag = false;
%         end
% 
% end % while; the flux convergence                                            
% 
% % Some debugging code
% loopTotal = loopTotal + loopCount;
% Max_Iterations = max(Max_Iterations,loopCount);
% Mean_Iterations_Per_Timestep = loopTotal/activeCount;
% if loopCount >= loopLimit;
%     disp('Max dQ exceeds tolerance: convergence not reached')
%     MaxdQ/Flux_Tolerance %#ok<NOPTS>
% end

%%%%%%%
            % Now that we know the final erosion rates, calculate a final elevation change based on the new timestep and dZdt.
            dZreg = Edot_reg .* dtReg;
            
                    
            if DETACH_LIMIT_STREAMS
               % Remove all mobile regolith wherever bedrock channels are defined
               dZreg(BRChan) = -Regolith_H(BRChan);
            end
            
            % Where has the regolith been completely removed?
            stripped = dZreg < -Regolith_H;
  
            dZreg(stripped) = -Regolith_H(stripped); % but only allow amount of regolith present to be eroded
  
            dzRegCum = dzRegCum + dZreg;    % Cumulative elevation change over external timestep
            topo = topo + dZreg;
            Regolith_H = Regolith_H + dZreg;
            Regolith_H(Regolith_H < 0) = 0;       % No negative thicknesses

            
            % Generate new regolith
            r_add = (rdot_H + rdot_dilate) .* dtReg; % generate this much regolith for the timestep

            topo = topo + rdot_dilate .* dtReg;      % add the dilation to the topography
            
            Regolith_H = Regolith_H + r_add;
            RegAdded = sum(sum(r_add));  % for tracking
            
            
            % This is the cumulative time of bedrock exposure to report 
            % to the bedrock erosion code
            exposed_time = exposed_time + (dtReg-dt_prime); 
            
%% Track mass conservation            
Reg_Vol = sum(sum(Regolith_H)) ;
Reg_Vol_init = Reg_Vol_init + RegAdded;
%Csed = SuspVol*rho_rock ./ AccGrid*dtReg; % For reference; kg sediment/m^3 water/yr


%% Sum flux into boundary cells that is removed from the grid
% Boundary_Loss = Boundary_Loss + dtReg * sum(sum( Qin(logical(BorderGrid)) )); 
% Boundary_Loss = (Boundary_Loss + dVperiodic);


end % internal while loop


