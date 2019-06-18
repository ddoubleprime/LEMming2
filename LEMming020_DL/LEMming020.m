% LEMming - a MATLAB-based landscape evolution model
% Copyright (C) 2007-2015 Dylan Ward
% 
% 
% Developer can be contacted at dylan.ward@uc.edu
% 
% and 
% 
% University of Cincinnati, 
% Dept. of Geology, 
% Cincinnati, OH 45221
% 
% 
% This program is free software; you can redistribute it and/or 
% modify it under the terms of the GNU General Public License as 
% published by the Free Software Foundation; either version 2 of 
% the License, or (at your option) any later version. See License.txt 
% and gpl-2.0.txt for the full terms of this license.
% -------------------------------------------------------------------

%% A MATLAB-based, heavily vectorized LEM with the ability to add layers of
%% variable lithology

%% D. Ward, 12/16/2007: Versions 002-015
%% D. Ward, 5/25/2011: Version 016
%% D. Ward, 4/29/2014: Version 020
%% V020 is a major upgrade, particularly in the regolith handling. Improvements include:
%%      - Significant improvements to regolith routing and transport
%%      - Bedrock channel erosion now works as expected

%% LEMming_020_DL - DETACHMENT-LIMITED version 5 Oct. 2015. 

% Variables should be specified in LEMming_Input_<setupname>.m
% or otherwise set by script before running LEMming



while t <= tmax       % MAIN TIME LOOP ---------------------------------%%

    % Get the state of precipitation for this timestep
        if ~UNIFORM_PRECIP
            % Rainmaker_lemLink
            TWI_lemLink
        end
        
    % Establish boundary conditions
        if OCEAN_BOUND
            OceanBoundary_lemLink;             % Handles variable sea-level BC
        else
            Enforce_Borders_lemLink;           
        end
    
    % See if rock uplift rate should change
        if VAR_ROCK_UPLIFT && ruindx < size(rock_uplift_mtx,1);
           if t >= rock_uplift_mtx(ruindx+1,1)
               ruindx = ruindx + 1;
               rock_uplift = rock_uplift_mtx(ruindx,2);
               U(:,:) = rock_uplift;       % Generic rock uplift rate field
           end
        end
        
    Streamgrid_lemLink;             % Update the flow matrices
    
         
    %%%% STRATIGRAPHY HANDLER %%%%
    
        % For each defined unit, Read top left and lower right x,y,z-values
        
            RTGrid(:,:) = 0;  % Reset rocktype grid - Default rocktype is zero
            depth_to_next(:,:) = inf; % same for depth grid - set inf depth unless specified below
            
            % Make an array of all unit IDs with a top bound above z_bound and
            % a bottom bound below max(topo), then use it, so the surface
            % intersection calculation doesn't need to go through the whole list
            % when there's lots of stratigraphy

            nearsurf = stratlist(7,:) <= max(max(topo)) & stratlist(4,:) >= z_bound;
            unitlist = stratlist(8,nearsurf);
      
        for unit = unitlist
        
            % RTGrid where topo(Xtl:Xlr,Ytl:Ylr) is less than top value and more than bottom value is assigned
            % the rock type index of that stratigraphic unit. If units are
            % defined to overlap, later units supercede earlier.
            rtype_here(:,:) = 0;
            topoabove = BRtopo(stratlist(3,unit):stratlist(6,unit),stratlist(2,unit):stratlist(5,unit)) <= stratlist(4,unit);
            topobelow = BRtopo(stratlist(3,unit):stratlist(6,unit),stratlist(2,unit):stratlist(5,unit)) >= stratlist(7,unit);                        
            rtype_here(stratlist(3,unit):stratlist(6,unit),stratlist(2,unit):stratlist(5,unit)) = topoabove .* topobelow;
            rtype_here = logical(rtype_here);
            
            RTGrid(rtype_here) = stratlist(1,unit);
            depth_to_next(rtype_here) = BRtopo(rtype_here) - stratlist(7,unit);

        
        end % for unit
        
        % Update the stratigraphy list z-values given the rock uplift rate
        %% Currently incompatible with spatially-varying rock uplift!
        stratlist(4,:) = stratlist(4,:) + (rock_uplift * dt);
        stratlist(7,:) = stratlist(7,:) + (rock_uplift * dt);
        
        % Now add in any surface layers tracked
        if DO_ROCKFALL
            RTGrid(RF_Debris_H > 0) = RF_Debris_Rtype;
            depth_to_next(RF_Debris_H > 0) = RF_Debris_H(RF_Debris_H > 0);
        end
        
    %%%% PROPERTY HANDLER %%%%
    
        % Establish grids for each property. Loop through units present in RTGrid 
        % and assign property values where RTGrid is each.       
       
        k(:,:) = k0;                % fluvial erodibility constant
        m(:,:) = m0;                % area exponent
        n(:,:) = n0;                % slope exponent
        kappa(:,:) = kappa0;        % regolith flux parameter, m/yr
        sc(:,:) = sc0;              % critical slope
        rfslope(:,:) = rfslope0;    % rockfall threshold slope
        
        for rtype = 1:numStypes;
        
            prtype_here = (RTGrid == rtype);

            k(prtype_here) = RTproplist(2,rtype);              % fluvial erodibility constant
            m(prtype_here) = RTproplist(3,rtype);              % area exponent
            n(prtype_here) = RTproplist(4,rtype);              % slope exponent
            kappa(prtype_here) = RTproplist(5,rtype);           % hillslope diffusivity, m/yr 
            sc(prtype_here) = RTproplist(6,rtype);              % critical slope for Roering nonlinear transport
            rfslope(prtype_here) = RTproplist(7,rtype);        % rockfall threshold slope
        
        end % for prtype
        
        
    % Need staggered grids for Kappa and sc. calculate mean values
    % at cell boundaries. edges are by default the kappa0 and sc0 values -
    % change to use nearest neighbor in relevant direction if problematic
    kappa_x(:,2:end-1) = ( kappa(:, 1:end-1) + kappa(:, 2:end) ) / 2;
    kappa_y(2:end-1,:) = ( kappa(1:end-1, :) + kappa(2:end, :) ) / 2;
    sc_x(:,2:end-1) = ( sc(:, 1:end-1) + sc(:, 2:end) ) / 2;
    sc_y(2:end-1,:) = ( sc(1:end-1, :) + sc(2:end, :) ) / 2;
    sThresh_x = sc_x-0.00001*sc_x;   % Prevents infinite erosion rates and backdiffusion at critical slope. 
    sThresh_y = sc_y-0.00001*sc_y;   % Prevents infinite erosion rates and backdiffusion at critical slope. 
    
    
    
    %%%% PROCESS HANDLER %%%%
            

            
        % Bedrock erosion
        switch BEDROCK_FLUV_MODE
            case 1
                Edot_bedrock = -k .* AccGrid.^m .* (sinSlopes).^n; % General stream-power model
            case 2
                Edot_bedrock = -k .* StreamTau.^n;    % Based on basal shear stress
            case 3
                Edot_bedrock = -k .* StreamPower;   % Based on total explicit stream power
            case 4
                Edot_bedrock = -k .* Omega;         % Based on unit stream power
            otherwise
                % Default to general stream-power model
                Edot_bedrock = -k .* AccGrid.^m .* (sinSlopes).^n;
        end
        

        % Hillslope erosion using a critical-slope (Roering) model
        % forward differences 
        dZdx(:,2:end-1) = ( topo(:,2:end) - topo(:,1:end-1) ) / dx ;
        dZdy(2:end-1,:) = ( topo(2:end,:) - topo(1:end-1,:) ) / dy ;
             
        % Apply slope limits to threshold slopes to prevent numerically infinite
        % transport and backdiffusion.
        dZdxL = dZdx;
        dzmask = (dZdx >= sThresh_x | dZdx <= -sThresh_x);
        dZdxL(dzmask) = sThresh_x(dzmask) .* sign(dZdx(dzmask));
        
        dZdyL = dZdy;
        dzmask = (dZdy >= sThresh_y | dZdy <= -sThresh_y);
        dZdyL(dzmask) = sThresh_y(dzmask) .* sign(dZdy(dzmask));
                         
        QregDiffX = (-kappa_x .* (dZdx))./(1-(abs(dZdxL)./sc_x).^2);    % Regolith transport flux, Roering 2001
        QregDiffY = (-kappa_y .* (dZdy))./(1-(abs(dZdyL)./sc_y).^2); 
        
        dQdx = ( QregDiffX( :, 2:end) - QregDiffX( :, 1:end-1) ) / dx;
        dQdy = ( QregDiffY( 2:end, :) - QregDiffY( 1:end-1, :) ) / dy;

        Edot_HS = (-dQdx -dQdy) ;
        Edot_HS(Edot_HS > 0) = 0;   % no deposition
        
        Edot_HS_fractional = Edot_HS .* (1-StreamFrac) ;  % correct for fraction of cell that is not channel

            % Correct for stream width and bedrock exposure
            Edot_bedrock(RNANS) = 0;
            Edot_bedrock_fractional = Edot_bedrock .* StreamFrac;
            
            % Calculate total erosion rate
            if EROSION_SUM_MODE == 1
                dHdtTot = Edot_HS_fractional + Edot_bedrock_fractional ; % Sum hillslope and channel components for each cell, proportionally to the fraction of the cell occupied by each
            elseif EROSION_SUM_MODE == 2
                dHdtTot = min(Edot_HS_fractional,Edot_bedrock_fractional) ; % Use maximum (negative) of hillsope or channel erosion, corrected for fraction of cell that is each
            else
                dHdtTot = min(Edot_HS,Edot_bedrock) ;   % Use maximum (negative), uncorrected for stream fraction of cell
            end
            
            Calc_dt_lemLink         % Calculate the dynamic timestep based on erosion rates.
            
            
        % Rockfall. Requires that dt be already calculated.
        if DO_ROCKFALL
            erodeMaxRate = RFerodeRate;
            RF20_lemLink;   % populates RFdz
        else
            RFdz(:,:) = 0;
        end
        
        dz = dt .* dHdtTot ;    % total elevation change this timestep due to erosion. 
        % RFdz is added to the topography by the landsliding routine
            
            
        % Rockfall layer magic. SHOULD rewrite to handle thickness changes
        % in RF code and just subtract otherwise eroded material here. At
        % least in DL version.
        if DO_ROCKFALL
            newDebris = RFdz > 0;
            eroded = dz + RFdz < 0;

            RF_Debris_H(newDebris) = RF_Debris_H(newDebris) + RFdz(newDebris);     % Add new rockfall debris into the debris thickness layer
            RF_Debris_H(eroded) = RF_Debris_H(eroded) + dz(eroded);   % Subtract eroded material
            RF_Debris_H(RF_Debris_H < 0) = 0;       % No negative thicknesses
            RF_Debris_H(BorderGrid) = 0;            % No storage on border

        end

                        
            % Update topography due to rock uplift and erosion
            topo = topo + dz + (U .* dt);
            
            
            % Enforce boundaries. Should recode enforce borders to go here.
            BRtopo = topo-RF_Debris_H;
            OERODE = BRtopo < z_bound;
                     
            topo(topo < z_bound) = z_bound;
            BRtopo(BRtopo < z_bound) = z_bound;

            
            
        % Track volumetric erosion rates by each process every 'tracktime' years
            if t_track >= tracktime || t == 0 || t >= tmax

                t_track = 0 + rem(t,tracktime);     % Reset tracking timer such that additional time due to 
                                                    % variable dt doesn't accumulate

                trackstep = trackstep + 1;

%                 TIMEtrackVol(trackstep) = t;
%                 FLUVtrackVol(trackstep) = (CellArea / dt) * sum(sum(dzFluvFrac));
%                 REGtrackVol(trackstep) = (CellArea / dt) * (sum(sum(RegAdded + sediment)));
%                 REGDZtrackVol(trackstep) = (CellArea / dt) * (sum(sum(dzRegCum + sediment)));
%                 RFtrackVol(1,trackstep) = (CellArea / dt) * sum(sum(RFdz(RFdz < 0)));
%                 RFtrackVol(2,trackstep) = (CellArea / dt) * sum(sum(RFdz(RFdz > 0)));
                
                % Check conservation of mass
                            
            
                % Keep track of total elevation change - should equal
                % init topo-final topo at the end of the run...
                dzCum = dzCum + dz;    

                
                Delta_Z = topo_init - (topo - rock_uplift .* t);   % Total bedrock elevation difference from initial condition
                Delta_V = sum(sum(Delta_Z)) * CellArea; % Volume of material that change represents
                
                Mass_Loss = sum(sum(dzCum)) * CellArea;  % Cumulative elevation change since the beginning of the run
                
                t %#ok<NOPTS>
                dt %#ok<NOPTS>

                Mass_Loss %#ok<NOPTS>
                

                % Volume of rock uplifted into the border area and removed
                RU_border = (rock_uplift * t * BorderArea);
                
%                ConservRatio = Boundary_Loss / ( dVS + dVN + dVE + dVW - RU_border) %#ok<NOPTS>
                
            end % if tracktime
        
        % Plot every "plottime" years
            if t_plot >= plottime || t == 0 || t >= tmax

                t_plot = 0 + rem(t,plottime);   % Reset plot timer such that additional time due to 
                                                % variable dt doesn't accumulate

                                
                
                stype = RTGrid;% .* (z_max_plot / numStypes);% + (topo/numStypes);
                coalcliffs = [0.3, 0.3, 0.3
                                1, 0.8, 0
                                0.3, 0.5, 0.3
                                0.3, 0.5, 0.3];

                figure(1); clf('reset'); 
                colormap(coalcliffs)
                surf(Xs,Ys,topo,stype+1,'FaceLighting','gouraud','FaceColor','interp',...
                   'AmbientStrength',ambient,'SpecularStrength',.2,'facealpha',1,'cdatamapping','direct'); 
                
                shading flat; brighten(0.5); lighting gouraud; light('Position',lightpos,'Style','infinite'); light('Position',[-10 -2 10],'Style','local','color',[1 .8 .6]); 
                set(gca,'DataAspectRatio',[1 1 1/VE]); view(az,el); zlim([z_bound z_max_plot]);
                title([run_name ' Topography and Lithology, Year ' int2str(t)]);
                


                qregmap = slopes;               
                qregmap(isinf(qregmap)) = max(max(qregmap(~isinf(qregmap))));
                figure(2); clf('reset'); 
                surf(Xs,Ys,topo,qregmap,'FaceLighting','gouraud','FaceColor','interp',...
                   'AmbientStrength',ambient,'SpecularStrength',.2,'facealpha',1,'cdatamapping','scaled');
                colorbar; colormap jet; hold on; caxis([min(min(qregmap)) max(max(qregmap))]); 
                set(gca,'DataAspectRatio',[1 1 1/VE]); view(az,el); zlim([z_bound z_max_plot]); 
                shading interp; 
                title([run_name ' slope, Year ' int2str(t)]);


%                 qregmap = StreamPower;               
%                 qregmap(isinf(qregmap)) = max(max(qregmap(~isinf(qregmap))));
%                 figure(3); clf; colormap jet; contour3(Xs,Ys,topo,cinterval:cinterval:z_max_plot,'k'); hold on; surfc(Xs,Ys,topo,qregmap); colorbar; set(gca,'DataAspectRatio',[1 1 1/VE]); view(-50,50); zlim([z_bound z_max_plot]); shading interp; 
%                 title([run_name ' Stream power, Year ' int2str(t)]);

                qregmap = t-RFAge;  
                if max(max(qregmap)) > min(min(qregmap))
                qregmap(isinf(qregmap)) = max(max(qregmap(~isinf(qregmap))));
                figure(4); clf('reset'); colormap jet; contour3(Xs,Ys,topo,cinterval:cinterval:z_max_plot,'k'); hold on; caxis([min(min(qregmap)) max(max(qregmap))]); surfc(Xs,Ys,topo,qregmap); colorbar; set(gca,'DataAspectRatio',[1 1 1/VE]); view(-50,50); zlim([z_bound z_max_plot]); shading interp; 
                title([run_name ' Rockfall debris age, Year ' int2str(t)]);
                end
                
%                 instedot = max((dz / dt),-.001);
%                 instedot = min(instedot,.001);
%                 figure(5); clf; colormap jet; contour3(Xs,Ys,topo,cinterval:cinterval:z_max_plot,'k'); hold on; caxis([min(min(instedot)) max(max(instedot))]); surfc(Xs,Ys,topo,instedot); colorbar; set(gca,'DataAspectRatio',[1 1 1/VE]); view(-50,50); zlim([z_bound z_max_plot]); shading interp; 
%                 title([run_name ' Instantaneous Erosion Rate (m/yr), Year ' int2str(t)]);
                
                drawnow

                if SAVEMODE == 1      % Save only the figure
                    disp('Saving current figure... ')

                    hgsave(['./' run_filename '/' run_name '_Plot' int2str(stateNo) '.fig'])

                    % Save the entire workspace the first time
                    if t == 0
                        save(['./' run_filename '/' run_name '_State' int2str(stateNo) '.mat']);
                    end

                    stateNo = stateNo + 1;

                    disp('Saved.')

                elseif SAVEMODE == 2    % Save the entire workspace
                    disp('Saving current state... ')

                    % Save the workspace
                    save(['./' run_filename '/' run_name '_State' int2str(stateNo) '.mat']);
                    stateNo = stateNo + 1;

                    disp('Saved.')

                end   % if SAVEMODE

                disp('Max erosion over last timestep: ')
                min(min(dz))
                disp('Max deposition over last timestep: ')
                max(max(dz))
                disp('Timestep: ')
                dt %#ok<NOPTS>
                disp('Time: ')
                t %#ok<NOPTS>
                
                if NoDistTargets ~= 0
                    disp([int2str(NoDistTargets) ' rockfall events had no distribution targets.'])
                    NoDistTargets = 0;
                end
                
            end % if t_plot

    %% Timestep is calculated by the regolith handler
        
        t = t + dt;             % Update time
        t_plot = t_plot + dt;   % Update plot time counter
        t_track = t_track + dt; % Update track time counter
        
end                % END MAIN TIME LOOP ---------------------------------%%

%% Call the finishing routine %%
FinishRun_lemLink


%%%%%%%% END OF LEMming.m %%%%%%%%