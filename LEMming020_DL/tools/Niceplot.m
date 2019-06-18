close all

load('/Users/dylanward/Desktop/LEMming Working/LEMming020_DL/RUNS/DL/slower 21-January-2018 01.08/slower_State1998.mat')

FigWidth = 800; % Specify the desired movie size
FigHeight = 600;    % in pixels
colordef white



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
z_max_plot = 100;
az = -50;
el = 50;

colormap(coalcliffs)

figure(1); clf; brighten(0.5);
surf(Xs,Ys,topo,stype+1,'FaceLighting','gouraud','FaceColor','interp',...
   'AmbientStrength',ambient,'SpecularStrength',.2,'facealpha',1,'cdatamapping','direct'); 

shading flat; lighting gouraud; light('Position',lightpos,'Style','infinite'); light('Position',[-10 -2 10],'Style','local','color',[1 .8 .6]); 
set(gca,'DataAspectRatio',[1 1 1/VE]); view(az,el); zlim([z_bound z_max_plot]); xlim([0 1200]);
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
%title([run_name ' Topography and Lithology, Year ' int2str(t)]);