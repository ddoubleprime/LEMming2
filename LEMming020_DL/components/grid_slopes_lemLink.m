% Figure out local slopes and flow directions

[R, slopes] = dem_flow(topo,dx,dy);

% RNANS = isnan(R);   % Find closed depressions for special treatment 
%%% CURRENTLY UNUSED: DEPRECATE?

% Negative slopes indicate closed depressions: set them to zero
% slope
slopes(slopes < 0) = 0;