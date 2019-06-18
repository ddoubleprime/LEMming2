% REGOLITH20_lemLink.m
% Component for LEMming to handle regolith production and transport by
% hillslope and fluvial processes


            % forward differences 
            dZdx(:,2:end-1) = ( topo(:,2:end) - topo(:,1:end-1) ) / dx ;
            dZdy(2:end-1,:) = ( topo(2:end,:) - topo(1:end-1,:) ) / dy ;
            
            
            BRtopof = fill_sinks(BRtopo);   % can't strip regolith below this limit
            
            
            %% Fix slope boundaries here if desired. ZERO BY DEFAULT


            %%%%%%%%% THIS IS THE REGOLITH HANDLER CODE %%%%%%%%%
            
        % Bedrock erosion
        switch BEDROCK_FLUV_MODE
            case 1
                Edot_bedrock = k .* AccGrid.^m .* (sinSlopes).^n; % General stream-power model
            case 2
                Edot_bedrock = k .* StreamTau.^n;    % Based on basal shear stress
            case 3
                Edot_bedrock = k .* StreamPower;  % Based on total explicit stream power
            otherwise
                % Default to general stream-power model
                Edot_bedrock = k .* AccGrid.^m .* (sinSlopes).^n;
        end
        
        % Correct for stream width and bedrock exposure
            Edot_bedrock = Edot_bedrock .* StreamFrac .* BRChan;
            

        % Regolith production rate
            rdot_H = rdot .* exp(-Regolith_H ./ rstar);     % Heimsathian
            rdot_H(Regolith_H <= Regolith_Prod_Thresh) = 0; % no regolith production where thinner than a threshold. Poor man's humped function.
            rdot_H(BorderGrid) = 0;
            
            
            % Scale regolith production by precipitation rate
            if SCALE_REG_PROD

                % Gradient of production rate fraction with rainfall rate
                rdotF = rdotF_P0 + (1-rdotF_P0) .* (PrecipGrid ./ RRateRef);
                rdot_H = rdot_H .* rdotF;

            end
            
            rdot_Tot = Edot_bedrock + rdot_H;
            rdot_dilate = ((rho_rock./rho_reg) - 1) .* rdot_Tot; % estimate of dilation rate of the regolith, used in timestep calculation. Actual dilation calculated based on production, below

            
        % Regolith transport
        
            % Active regolith thickness is difference between topography and
            % pit-filled bedrock top. Needs to deal with rockfall layers
            % here.
            Regolith_H_active = topo - BRtopof + RF_Debris_H;
            
            % Enforce mobile regolith thickness
            Mobile_H = min(Regolith_H_active,Max_Mobile_H);
            Mobile_H(BorderGrid) = Max_Mobile_H;
            
            % Mobile thickness at cell edges
            MobH_X_int = (Mobile_H(:,1:end-1) + Mobile_H(:,2:end)) / 2;  % interior grid, 0x1 smaller than topo
            MobH_X(:,2:end-1) = MobH_X_int;

            MobH_Y_int = (Mobile_H(1:end-1,:) + Mobile_H(2:end,:)) / 2;  % interior grid, 1x0 smaller than topo
            MobH_Y(2:end-1,:) = MobH_Y_int;
            

        Qout = kt .* AccGrid.^mt .* sign(slopes) .* (sinSlopes).^nt .* StreamWs;    % Regolith transport flux by stream transport
        Qout(Regolith_H <= 0) = 0;      % no sediment, no sediment discharge
        Qout(BorderGrid) = 0;         % Trap all regolith on boundaries for accounting

        % Sum local fluxes. Relies on LocalQ.mex being compiled and in the
        % path.
        Qin = LocalQ(x,y,Qout,Ts);

        % Sum contributions to the flux divergence
        Edot_Fluv = (Qin - Qout) ./ (dx*dy);
        
        
        % Apply slope limits to threshold slopes to prevent numerically infinite
        % transport.
        dZdxL = dZdx;
        dZdxL(dZdx >= sThresh | dZdx <= -sThresh) = sThresh;
        dZdyL = dZdy;
        dZdyL(dZdy >= sThresh | dZdy <= -sThresh) = sThresh;
                        
        UregDiffX = (-kappa .* (dZdx))./(1-(abs(dZdxL)./sc).^2);    % Regolith transport rate, Roering 2001
        UregDiffY = (-kappa .* (dZdy))./(1-(abs(dZdyL)./sc).^2); 
        
        overflux = abs(UregDiffX * dt) > dx;
        UregDiffX(overflux) = sign(UregDiffX(overflux)) .* (dx/dt);     % limit flux rate to traversing one cell
        overflux = abs(UregDiffY * dt) > dy;
        UregDiffY(overflux) = sign(UregDiffY(overflux)) .* (dy/dt);     % limit flux rate to traversing one cell
        
        % Sum contributions to regolith motion
            
            UregX = UregDiffX; % + any other processes you might add: stochastic mass wasting, superdiffusive flux, etc.
            UregY = UregDiffY;


 
        % Enforce internal boundaries to ensure that no regolith is
        % transported from bare rock. From GC2D...
            CLASS = ( Regolith_H > 0 ) ;
    
            DCLASSx = zeros(y,x+1);
            DCLASSy = zeros(y+1,x);
            
            DCLASSx(:,2:end-1) = ( CLASS(:,2:end) - CLASS(:,1:end-1) ) .* ...
                        sign( dZdx(:,2:end-1) ) ;
                     
            DCLASSy(2:end-1,:) = ( CLASS(2:end,:) - CLASS(1:end-1,:) ) .* ...
                        sign( dZdy(2:end-1,:) ) ;
            
            indCx = ( DCLASSx == -1 ) ;
            UregX(indCx) = 0;
                        
            indCy = ( DCLASSy == -1 ) ;
            UregY(indCy) = 0;
            
  
            
        % Diffusive regolith flux
            QregX = (UregX .* MobH_X);  
            QregY = (UregY .* MobH_Y);
            

        % Convert to erosion rates
            dQdx = ( QregX( :, 2:end) - QregX( :, 1:end-1) ) / dx;
            dQdy = ( QregY( 2:end, :) - QregY( 1:end-1, :) ) / dy;

            Edot_reg = (-dQdx -dQdy) .* (1-StreamFrac) ;

            
                      % Figure out a timestep. Don't let slopes reverse. %
           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                   
        %%% CALCULATE TIMESTEP
        
        if ( VARIABLE_DT_TOGGLE )
            
            %% now that we know the rate of change in surface heights due to  
            %% regolith motion and production we need to know over 
            %% what period of time we can project forward with these rates and 
            %% maintain stability of the land surface. The basic idea here is that
            %% we don't want to take a timestep any longer then it would take to 
            %% reverse the land surface slope between two cells, such that regolith 
            %% should be flowing in the other direction. In fact, let's make our 
            %% timestep much less than that.
            
            %% this calculation sets the timestep such that the change
            %% in surface elevation nowhere exceeds a set fraction
            %% of the local standard deviation in surface elevations
            
            
            dHdtTot = Edot_reg + Edot_Fluv + rdot_dilate ;
            adHdt = abs(dHdtTot) ;
            
            adZdt = abs(Edot_bedrock) ;
            
            % something like standard deviation of 3x3 cell areas around each cell
            topoPad(2:end-1,2:end-1) = topo;
            topoMean = filter2( dtfilt, topoPad, 'valid' ) ;
            dHmax = sqrt( filter2( dtfilt, (topoMean - topo).^2 ) ) ;
            
            topoPad(2:end-1,2:end-1) = BRtopo;
            topoMeanbr = filter2( dtfilt, topoPad, 'valid' ) ;
            
            dHmaxbr = sqrt( filter2( dtfilt, (topoMeanbr - BRtopo).^2 ) ) ;
            
            % find limiting timestep for each considered cell
            ind = ( adHdt~=0 & dHmax~=0 & ~dtMask) ;    
            dtLimits = dHmax(ind)./adHdt(ind) ;
            [dtH, idt] = min( dtLimits ) ;
            
            ind = ( adZdt~=0 & dHmax~=0 & ~dtMask) ;    
            dtLimits = dHmaxbr(ind)./adZdt(ind) ;
            [dtZ, idtz] = min( dtLimits ) ;
            
            
            % Separately find the limiting timestep for regolith production
            % accuracy.             
            % Regolith thickness increase should not exceed a preset amount in a
            % timestep or the production function calculation becomes
            % inaccurate. - DEPRECATE 7/9/15 - exact regolith calculation
            % implemented
            
            % dtR = min(min(Reg_prod_limit ./ ( rdot_Tot(rdot_Tot >= 0) + rdot_dilate(rdot_Tot >= 0) ) ));
            
            % limit timestep to dtR or some fraction of the calculated timestep
            %if ~isempty(dtR)
            %    dt = min(dtH/4,dtR);
           % else
                dt = dtH/4;
           % end
            if ~isempty(dtZ)
                dt = min(dtZ/2,dt);
            end
            dt = min(dt,dtMax);
            dt = max(dt,dtMin);
            % catch an error
            if isempty(dt)  
                dt = dtDefault ;
            end
            
        else

            dt = dtDefault ;
                        
        end %- VARIABLE_DT_TOGGLE -%
            
        
        
            SuspLoss = FineFraction .* Regolith_H_active .* StreamFrac;  % This is crude and does not account for coarsening of the regolith in streams, so use with caution.
            
            % Analytical solution for exponential regolith production over
            % the timestep, adding in bedrock erosion contribution
           
            H_prod = rstar .* log(((rdot_H .* dt) ./ rstar) + 1) + Edot_bedrock * dt;
            H_dilate = ((rho_rock./rho_reg) - 1) .* (H_prod + Edot_bedrock * dt);
            
            dZreg = (Edot_reg + Edot_Fluv) .* dt + H_dilate - SuspLoss;
            dHreg = (Edot_reg + Edot_Fluv) .* dt + H_prod + H_dilate - SuspLoss;
                        
            % subtract out reg production in channel portion of cell
            dHreg(StreamStripped) = dHreg(StreamStripped) .* (1 - StreamFrac(StreamStripped));
            
            % find cells completely drained of regolith
            stripped = dHreg < -Regolith_H_active;
            
            % Track inaccuracy due to pulling regolith from cells stripped mid-timestep. Could be used to adjust the timestep.
            Strip_Error = Strip_Error + sum( dHreg(stripped) - Regolith_H_active(stripped) ); 
 
            % Correct elevation/thickness change so as to avoid bedrock
            % erosion by regolith transport
            dZreg(stripped) = -Regolith_H_active(stripped);
            dHreg(stripped) = -Regolith_H_active(stripped);
            
            % Here figure out cells in which the fluvial erosion exceeded
            % the stream fraction of the cell's total regolith, even if the cell
            % was not fully stripped. We'll treat these as having bedrock
            % channels.
            dHFluv = (Edot_Fluv .* dt) - SuspLoss;
            StreamStripped = ( dHFluv < StreamFrac .* -Regolith_H );       
            StreamStripped = StreamStripped | stripped;         % if the cells were stripped fully, include them.
            
            if ~DETACH_LIMIT_STREAMS
                BRChan = StreamStripped | Regolith_H < Max_Mobile_H;  % if regolith is all mobile within the stream fraction of the cell, it should be a detachment-limited channel, yes?
            end
            
            % Update topography and regolith thicknesses
            topo = topo + dZreg;
            Regolith_H = Regolith_H + dHreg;
            
            % Now, avalanche off of slopes greater than the angle of repose
            % for regolith
            skipRFflag = 0;  % set a flag so the rockfall routine doesn't fire on the same timestep. may artifically limit rockfall rate at small values of avalanchFreq.
        if ( AVALANCH_TOGGLE && ( rand < dt*avalanchFreq ) )

            LSMODE = 0;     % zero for regolith, one for rockfall       
            ShallowLS_lemLink
            skipRFflag = 0; % zero here enbales rockfall regardless
            
        end

