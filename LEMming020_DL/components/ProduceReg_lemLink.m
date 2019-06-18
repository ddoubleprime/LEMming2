%% ProduceReg_lemLink.m - LEMming component that figures out new regolith
%% production rates for each timestep, based on current thicknesses and
%% precipitation (infiltration, eventually) rates

%% v.016 for LEMming v.016: Rolled into REGOLITH20 for v. 020

% Rates according to exponential decline scheme
rdot_H = rdot .* exp(-Regolith_H ./ rstar);     % Heimsathian
%rdot_H = rdot_H .* (1-(StreamWs ./ GridDelta));    % Don't produce where there are streams
rdot_H(BRChan) = 0; % don't produce regolith in forced-bedrock channels
%rdot_H(BorderGrid) = 0;

% Scale regolith production by precipitation rate
if SCALE_REG_PROD
   
    % Gradient of production rate fraction with rainfall rate
    rdotF = rdotF_P0 + (1-rdotF_P0) .* (PrecipGrid ./ RRateRef);
    rdot_H = rdot_H .* rdotF;

end

rdot_dilate = ((rho_rock./rho_reg) - 1) .* rdot_H; % dilation rate of the regolith