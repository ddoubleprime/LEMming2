% LEMming_CliffTrax.m - extracts 

clc
clear all
close all

filename = 'LessSteepToFail'; % base filename
process_folder = '/Users/dylanward/Desktop/LEMming Working/LEMming020_DL/RUNS/DL/LessSteepToFail 11-April-2018 09.53';

varname = 'topo'; % variable to extract along with cliff retreat rate

%FigWidth = 700; % Specify the desired movie size
%FigHeight = 700;    % in pixels

MAX_FLAG = 1;   % Extract maximim of varname in each swath
MIN_FLAG = 1;   % Extract minimum of varname in each swath
MEAN_FLAG = 1;  % Extract maximim of varname in each swath

MASK_FLAG = 1;  % mask the grid for calculation by a threshold value
mask_val = 0;

% Color axis
clim0 = 0; 
clim1 = 1e2;

cmap = 'jet';

% swath subdomain (cells) - set above image size for whole area
x0 = 40;
x1 = 80;
y0 = 1;%+150;
y1 = 2000;




%%% End of inputs %%%

cd '/Users/dylanward/Desktop/LEMming Working/LEMming020_DL'

filelist = dir([process_folder '/' filename '_State*.mat']);

[numfiles, a] = size(filelist);

s = cell(numfiles,1);

for fn = 1:numfiles
   
   currfn = [process_folder '/' filename '_State' int2str(fn-1) '.mat'];
   s{fn} = struct(load(currfn,varname,'t','topo','rock_uplift','dy'));
   
end


    [yimg, ximg] = size(s{2}.(varname));
    y1 = min(y1,yimg);
    x1 = min(x1,ximg);

    % allocate results array. 
    % Columns: 1,t; 2,rock_uplift; 3,retreat rate; 4,min(varname); 5,mean(varname); 6,max(varname)
    results = -9999 * ones(numfiles,6);
    cliffpos = zeros(numfiles,1);
    for frame = 1:numfiles %#ok<*UNRCH>
        
        % store read-in scalar variables
        results(frame,1) = s{frame}.t;
        results(frame,2) = s{frame}.rock_uplift;
        
        % copy gridded variables for manipulation
        img = s{frame}.(varname);
        topo = s{frame}.topo;
        
        % calculate cliff backwearing rate
        topotrim = topo(y0:y1,x0:x1);
        topoavg = mean(topotrim,2);
        cliffpos(frame) = mean(find(topoavg >= 0.9*max(topoavg)) * s{frame}.dy);
        
        results(frame,3) = cliffpos(frame);
        
        % other variables of interest
        
        if strcmp(varname,'RFAge')
            img = s{frame}.t-img;
        end
        
        imgtrim = img(y0:y1,x0:x1);
                
        %plot for reference
%         figure(1); set(gcf,'Position',[1 1 FigWidth+1 FigHeight+1]);

%         imagesc(imgtrim);
%         axis image; colormap(cmap); colorbar; caxis([clim0 clim1]); title(['Year ' int2str(s{frame}.t)]);
%         figure(1)
%         pause(.02)
        
        if MASK_FLAG
            imgmask = imgtrim(abs(imgtrim) > mask_val);
        else
            imgmask = imgtrim;
        end
        
        if ~isempty(imgmask)
            
            if MIN_FLAG
                minimg = min(min(imgmask));
                results(frame,4) = minimg;
            end

            if MEAN_FLAG
                meanimg = mean(mean(imgmask));
                results(frame,5) = meanimg;
            end

            if MAX_FLAG
                maximg = max(max(imgmask));
                results(frame,6) = maximg;
            end
            
        end
        
    end
    
figure(2); clf
hold on
subplot(3,1,1)
plot(results(:,1),results(:,2),'k')
ylabel('uplift rate m/y')
subplot(3,1,2)
hold on
presults=results; % plotting temp
presults(results == -9999) = nan;
plot(presults(:,1),presults(:,6),'r')       % max
plot(presults(:,1),presults(:,5),'g')       % mean
plot(presults(:,1),presults(:,4),'b')       % min
subplot(3,1,3)
plot(presults(:,1),cliffpos,'b')    % cliff position in y direction
xlabel('time, yr')
ylabel('cliff position, y coordinate (m)')

csvwrite(['results_' filename '_' varname '.csv'],results)


