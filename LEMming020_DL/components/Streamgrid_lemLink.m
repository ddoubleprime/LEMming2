% Streamgrid_lemLink.m - LEMming component that recalculates the Dinfinity
% flow grid and other stream network parameters as needed by each process.

%% DETACHMENT-LIMITED v020

          topof = fill_sinks(topo);  

        % Figure out local slopes and flow directions

            [R, slopes] = dem_flow(topof,dx,dy);


        % Negative slopes indicate closed depressions: set them to zero
        % slope
            slopes(slopes < 0) = 0;

            RNANS = isnan(R);   % Find closed depressions for special treatment
            slopes(RNANS) = 0;
            slopes(BorderGrid) = 0;
            sinSlopes = sin(atan(abs(slopes)));
            

        % Accumulation areas
            Ts = flow_matrix(topof,R);
            AccGrid = upslope_weight(topof,Ts,PrecipGrid);
            AccGrid = (dx * dy) * postprocess_plateaus(AccGrid,topof);

        % Darcy-Weisbach flow velocity
            
            % Uses dimensionless slope. No zeros allowed in D-W calculation
            Wslopes = slopes;
            Wslopes(Wslopes == 0) = min(min(Wslopes(Wslopes > 0)));
        
            Qw = (AccGrid/spyr);    % Runoff in m^3/sec
            StreamUs = ((8*g * sqrt(stream_WDR * Qw) .* Wslopes ) ./ (stream_WDR * f)).^(2/5);    
            
        % Stream widths
            StreamWs = sqrt( (stream_WDR .* Qw) ./ StreamUs);

        % Fraction of each cell occupied by a stream
            StreamFrac = StreamWs / GridDelta;

        % Stream Depths
            StreamDs = Qw ./ (StreamUs .* StreamWs);

        % Shear Stress
            % NEEDS SIN SLOPES for steep slopes
            StreamTau = rho_w * g .* StreamDs .* sinSlopes;            

        % Raw stream power
            StreamPower = StreamTau .* StreamUs;
            
        % Unit stream power
            Omega = StreamPower ./ StreamWs;
        
        % Correct stream widths after depth/stress calculation
            StreamWs(StreamWs > dx) = dx; % no channel should
                                          % be wider than a grid cell for
                                          % all other calculations. If you
                                          % have a lot of channels wider
                                          % than a grid cell, consider
                                          % running with a larger grid spacing.
