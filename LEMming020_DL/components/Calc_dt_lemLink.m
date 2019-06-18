%% Calc_dt_lemLink.m - Calculates the dynamic timestep for the detchment-limited version of LEMming.


            
                      % Figure out a timestep. Don't let slopes reverse. %
           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                   
        %%% CALCULATE TIMESTEP
        
        if ( VARIABLE_DT_TOGGLE )
            
            %% now that we know the rate of change in surface heights due to  
            %% erosion we need to know over 
            %% what period of time we can project forward with these rates and 
            %% maintain stability of the land surface. The basic idea here is that
            %% we don't want to take a timestep any longer then it would take to 
            %% reverse the land surface slope between two cells, such that regolith 
            %% should be flowing in the other direction. In fact, let's make our 
            %% timestep much less than that.
            
            %% this calculation sets the timestep such that the change
            %% in surface elevation nowhere exceeds a set fraction
            %% of the local standard deviation in surface elevations
            
            
            adHdt = abs(dHdtTot) ;
     
            
            % something like standard deviation of 3x3 cell areas around each cell
            topoPad(2:end-1,2:end-1) = topo;
            topoMean = filter2( dtfilt, topoPad, 'valid' ) ;
            dHmax = sqrt( filter2( dtfilt, (topoMean - topo).^2 ) ) ;
            
            
            % find limiting timestep for each considered cell
            ind = ( adHdt~=0 & dHmax~=0 & ~dtMask ) ;    
            dtLimits = dHmax(ind)./adHdt(ind) ;
            [dtH, idt] = min( dtLimits ) ;
                        
            dt = dtH/dt_safety_factor;


            dt = min(dt,dtMax);
            dt = max(dt,dtMin);
            % catch an error
            if isempty(dt)  
                dt = dtDefault ;
            end
            
        else

            dt = dtDefault ;
                        
        end %- VARIABLE_DT_TOGGLE -%
