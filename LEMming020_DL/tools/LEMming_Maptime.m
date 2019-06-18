% LEMming_Maptime.m - a different take on the movie maker. Load all LEM
% outputs from a folder, plot the variable of interest sequentially as a
% movie in map view.

clc
clear all
close all

filename = 'lesstilt'; % base filename
process_folder = '/Volumes/MacintoshHD/Users/dylanward/Desktop/Dropbox/Work/matlab_shared/LEMming020_DL/RUNS/DL/lesstilt 10-October-2015 20.51';

varname = 'Edot_HS';

FigWidth = 700; % Specify the desired movie size
FigHeight = 700;    % in pixels

LOG_FLAG = 0;   % Plot log10(varname). Good for things like discharge
DIFFS_FLAG = 0; % Plot differences between sucessive outputs rather than values
MASK_FLAG = 0;  % Plot location of values above a threshold (can be used with values or diffs)

mask_val = 1e3;

% Color axis
clim0 = -2e-4; 
clim1 = 0;

cmap = 'jet';

% plot subdomain (cells) - set above image size for whole area
x0 = 1;%+200;
x1 = 1000;
y0 = 1;%+150;
y1 = 2000;




%%% End of inputs %%%

cd /Volumes/MacintoshHD/Users/dylanward/Desktop/Dropbox/Work/matlab_shared/LEMming020_DL/

filelist = dir([process_folder '/' filename '_State*.mat']);

[numfiles, a] = size(filelist);

s = cell(numfiles);

for fn = 1:numfiles
   
   currfn = [process_folder '/' filename '_State' int2str(fn-1) '.mat'];
   s{fn} = struct(load(currfn,varname,'t'));
   
end
        [yimg, ximg] = size(s{2}.(varname));
        y1 = min(y1,yimg);
        x1 = min(x1,ximg);

if DIFFS_FLAG

    for frame = 2:numfiles

       figure(1); set(gcf,'Position',[1 1 FigWidth+1 FigHeight+1]);

        if LOG_FLAG
            img = log10(s{frame}.(varname))-log10(s{frame-1}.(varname));
        else
            img = s{frame}.(varname)-s{frame-1}.(varname);
        end

        if MASK_FLAG
            img = img >= mask_val;
        end
            
        imgtrim = img(y0:y1,x0:x1);
        imagesc(imgtrim);
        axis image; colormap(cmap); colorbar; caxis([clim0 clim1]); title(['Year ' int2str(s{frame-1}.t) ' - ' int2str(s{frame}.t)]);
        figure(1)
        pause(1)
        M(frame-1) = getframe(gcf);  %#ok<SAGROW> % Capture movie frame
        mname = ['_Delta' varname];
    end
    
else  % ~DIFFS_FLAG
    
    for frame = 1:numfiles %#ok<*UNRCH>

       figure(1); set(gcf,'Position',[1 1 FigWidth+1 FigHeight+1]);

        if LOG_FLAG
            img = log10(s{frame}.(varname));
        else
            img = s{frame}.(varname);
        end
        
        if MASK_FLAG
            img = img >= mask_val;
        end
        
        imgtrim = img(y0:y1,x0:x1);
        imagesc(imgtrim);
        axis image; colormap(cmap); colorbar; caxis([clim0 clim1]); title(['Year ' int2str(s{frame}.t)]);
        figure(1)
        pause(1)
        M(frame) = getframe(gcf);  %#ok<SAGROW> % Capture movie frame
        mname = varname;
        
    end
    
end
% movie(gcf,M)

movie2avi(M,[process_folder '/' filename '_' mname])