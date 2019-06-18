% MakeLEMmingMov.m - restores each saved step in the LEMming model run, plots it,
% captures the frame, and makes an .avi movie from it.

clc
close all

FigWidth = 800; % Specify the desired movie size
FigHeight = 600;    % in pixels
startstate = 0;     % 0 is default. N-1 starting state to plot (if only want to plot from middle of run)
nth_frame = 1;      % integer, plot every nth_frame
colordef none

if exist('M','var'); clear M; end  % Reinitialize the movie matrix

if ~exist('run_name','var')
    run_name = ' ';      % Copy the run name here or open any statefile from the run before running MakeLEMmingMov.m
    run_filename = ' ';  % Copy the run filename (Folder name) here or "
end


% Load the initial state 
try     % in case something isn't right in the run folder
    
    % Load the final workspace and read the SAVEMODE and final
    % stateNo
    %cd ~/Desktop/Dropbox/Work/MATLAB_SHARED/LEMming/LEMming020_DL
    load(['./' run_filename '/' run_name '_EndState.mat'],'SAVEMODE','stateNo');

    stateNo = stateNo - 1; % Because it would have been incremented past the final state
    
    frame = 1;  % initial frame
    
    % Loop through states and load each file
    for state = startstate:nth_frame:stateNo
                                
            if SAVEMODE == 1      % Saved only the figures

                close all
                open(['./' run_filename '/' run_name '_Plot' int2str(state) '.fig'])
                set(gcf,'Position',[1 1 FigWidth+1 FigHeight+1])
                           
            elseif SAVEMODE == 2    % Saved the entire workspaces


                % Load the workspace
                load(['./' run_filename '/' run_name '_State' int2str(state) '.mat']);

                % run_filename override - if files were mixed from
                % different runs
                %run_filename = '/Volumes/MacintoshHD/Users/dylanward/Desktop/Dropbox/Work/matlab_shared/LEMming020_dev/RUNS/DEV/50mrockfall';
                
                
                % Remake the plot
               stype = (RTGrid .* (z_max_plot / numStypes)) + (topo/numStypes);
               figure(1); clf; set(gcf,'Position',[1 1 FigWidth+1 FigHeight+1]); 
                               stype = RTGrid;% .* (z_max_plot / numStypes);% + (topo/numStypes);
                coalcliffs = [0.3, 0.3, 0.3
                                1, 0.8, 0
                                0.3, 0.5, 0.3
                                0.3, 0.5, 0.3];
                            
                              % shaded relief plot controls
                lightpos = [10 5 10];
                ambient = 0.5;
                VE2 = 2;
                z_max_plot = 300;
                az = -50;
                el = 50;

                colormap(coalcliffs)
                
                figure(1); clf; brighten(0.5);
                surf(Xs,Ys,topo,stype+1,'FaceLighting','gouraud','FaceColor','interp',...
                   'AmbientStrength',ambient,'SpecularStrength',.2,'facealpha',1,'cdatamapping','direct'); 
               
                shading flat; lighting gouraud; light('Position',lightpos,'Style','infinite'); light('Position',[-10 -2 10],'Style','local','color',[1 .8 .6]); 
                set(gca,'DataAspectRatio',[1 1 1/VE]); view(az,el); zlim([z_bound z_max_plot]);
                title([run_name ' Topography and Lithology, Year ' int2str(t)]);
               
%                colormap jet; contour3(Xs,Ys,topo,cinterval:cinterval:z_max_plot,'k'); hold on; surfc(Xs,Ys,topo,stype); set(gca,'DataAspectRatio',[1 1 1/VE]); view(-50,50); zlim([z_bound z_max_plot]); shading interp;
%                title([run_name ' Topography and Lithology, Year ' int2str(t)]);

%                 regmap = Regolith_H; 
%                 RegHCutoff = min( max(max(Regolith_H)), 3 );   % Upper color scaling bound
%                 minreg = min(min(Regolith_H));  % lower color scaling bound
%                 % Trap for a plotting crash that happens if RegHCutoff and
%                 % minreg are the same
%                 if RegHCutoff == minreg; RegHCutoff = minreg + 0.1*minreg; end
%                 
%                 % override
%                 %minreg = 0; RegHCutoff = 1;
%                 %z_max_plot = 100;
%                 
% %                 figure(1); clf; set(gcf,'Position',[1 1 FigWidth+1 FigHeight+1]); colormap jet; contour3(Xs,Ys,topo,cinterval:cinterval:z_max_plot,'k'); hold on; caxis([minreg RegHCutoff]); surfc(Xs,Ys,topo,regmap); set(gca,'DataAspectRatio',[1 1 1/VE]); view(-50,50); zlim([z_bound z_max_plot]); shading interp;
% %                  title([run_name ' Regolith Thickness (m), Year ' int2str(t)]); colorbar;
                 
                drawnow
                % ADDITIONAL PLOTS we may wish to make
                %             figure(2); imagesc(dz); view(2); shading flat; colorbar;
                %             figure(3); imagesc(slopes); view(2); shading flat; colorbar;
                %             figure(4); imagesc(curves); view(2); shading flat; colorbar;
                %             figure(5); imagesc(AccGrid); view(2); shading flat; colorbar;
                %             drawnow
                
            end   % if SAVEMODE

        % Get the movie frame from the state. If multiple plots, this code
        % needs to be put inside the IF SAVEMODE bit after each plot
        % command and separate M variables assigned.
        pause(0.1); % Short pause to make sure everything draws before the frame is captured
        M(frame) = getframe(gcf);  %#ok<SAGROW>
        frame = frame + 1;

    end % for state

    disp 'Saving movie...'
    MOVERROR = 0;
    
    mov_name = run_name;
    if isempty(mov_name); mov_name = 'LEMming_Movie'; end
    
    try
        vfile = VideoWriter(['./' run_filename '/' mov_name '.avi']);
        vfile.FrameRate = 10;
        open(vfile)
        writeVideo(vfile,M)
        close(vfile)
        %movie2avi(M,['./' run_filename '/' mov_name '.avi'],'compression','none','fps',10)
    catch movermsg
        MOVERROR = 1;
    end

    if ~MOVERROR
        % Prompt to toggle deletion of intermediate files to free disk space. Does not
        % remove Year0 or EndState.
        DO_CLEANUP = input('Remove all but the initial and final state savefiles? Y/N [Y]: ', 's');
        if isempty(DO_CLEANUP) || DO_CLEANUP == 'y';
            DO_CLEANUP = 'Y';
        end

        if DO_CLEANUP == 'Y'
           for state = 1:stateNo 
               if SAVEMODE == 1
                    delete(['./' run_filename '/' run_name '_Plot' int2str(state) '.fig'])
               elseif SAVEMODE ==2
                    delete(['./' run_filename '/' run_name '_State' int2str(state) '.mat'])
               end
           end

        end
    end
    disp 'Done!'
catch
        disp('One or more required files missing or run name not specified. Movie not made.')
end

if MOVERROR
    disp 'Error converting movie to AVI. Movie not saved.'
    disp(movermsg)
end

% END of MakeLEMmingMovie.m