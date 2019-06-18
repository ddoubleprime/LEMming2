
        
        % Add in suspended transport --  NOT FUNCTIONAL
        % Proportion of suspended sed scales with discharge
        % UregSusp = SuspRef .* (UregFluv ./ UregSuspRef); 

%         % Suspended flux is constant proportion of bedload flux
%         QregSusp = SuspRef .* QregFluv; 


% %        Settling velocity after Ferguson and Church (2006), assuming seive
% %        diameters and standard values for gravity and viscosity
%         Ws = 16.17*S50^2 / (1.8e-5 + sqrt(12.1275*S50^3));
%         
%         Rouse = Ws./(0.4 * sqrt(StreamTau./rho_w));  % Rouse number
%         % Rouse numbers below 2.5 indicate partial suspension. So, let's
%         % simply suspend the fine fraction in proportion to the amount that
%         % the Rouse number is below 2.5.
%         SuspFrac = max(0,1 - (Rouse / 2.5)) ;
%         SuspVol = SuspFrac .* Mobile_H .* FineFraction .* StreamWs;
%         Mobile_H = Mobile_H - (SuspVol/CellArea);   % Remove suspended sediment from the mobile layer
%         
%         % Rouse numbers below 0.8 indicate wash load; this can simply be
%         % removed from the grid and no longer tracked, loosening stability
%         % requirements.
%         WashVol = max(0,1 - (Rouse / 0.8)) .* SuspVol;
%         SuspVol = SuspVol - WashVol;    % Remove wash load from the tracked transport flux
%         WashQ = WashVol .* (StreamUs*spyr);
        
        % QregSusp = SuspVol .* (StreamUs*spyr); % Cubic meters suspended load flux per year 