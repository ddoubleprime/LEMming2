% ShallowLS_lemLink.m      - LEMming component that uses a stochastic downhill
% avalanching routine (from Kessler's GC2D avalanching) to handle rapid nonlocal
% fluxes of regolith. Part of Regolith20.




        
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% AVALANCHE REGOLITH OFF OF STEEP SURFACES
        

            
            %% move regolith downslope until the land surface is everywhere
            %% less than or near the angle of repose
            
            if LSMODE == 1;             % rockfall
                Ho = max( RFH, 0 ) ;
                dHRepose = dHReposeRF;
            else                        % regolith
                Ho = max( Regolith_H, 0 ) ;
                dHRepose = dHReposeReg;
            end
            
            avctMax=1e4;
            for avct = 1:avctMax

                % elevation differences, west-east                
                dZidx_down(:,2:end) = max( 0, ( topo(:,2:end) - topo(:,1:end-1) ) ) ;
                % elevation differences, east-west
                dZidx_up(:,1:end-1) = max( 0, ( topo(:,1:end-1) - topo(:,2:end) ) ) ;
                % take positive elevation differences
                dZidx = max( dZidx_up, dZidx_down ) ;

                % elevation differences, south-north
                dZidy_left(2:end,:) = max( 0, ( topo(2:end,:) - topo(1:end-1,:) ) ) ;
                % elevation differences, north-south
                dZidy_right(1:end-1,:) = max( 0, ( topo(1:end-1,:) - topo(2:end,:) ) ) ;
                % take positive elevation differences
                dZidy = max( dZidy_left, dZidy_right ) ;

                % elevation difference magnitude
                grad = sqrt( dZidx.^2 + dZidy.^2 );
                % total elev diffs toward cells
                gradT = dZidy_left + dZidy_right + dZidx_down + dZidx_up ;
                % set zero to 1 for division
                gradT(gradT==0) = 1;
                % ignore steep slopes with very thin mobile material
                grad(Ho < 0.001) = 0 ;

                mxGrad = max(max( grad ) ) ;
        
                if ( mxGrad <= 1.1*dHRepose )
                    break ;
                end

                % Amount of adjustment - 1/3 elevation difference minus
                % height equivalent of AOR at current cellsize (dHRepose)
                delH = max( 0, ( grad - dHRepose ) / 3 ) ;
        
                % subtract adjustment from deposit thickness and enforce positive-definite
                Htmp = Ho ;        
                Ho = max( 0, Htmp - delH );
                % total final change in thickness
                delH = Htmp - Ho ;
        
                % change in thickness due to west->east motion
                delHup(:,2:end) = delH(:,1:end-1) .* dZidx_up(:,1:end-1)./gradT(:,1:end-1) ;
                % change in thickness due to east->west motion 
                delHdn(:,1:end-1) = delH(:,2:end) .* dZidx_down(:,2:end)./gradT(:,2:end) ;
                % change in thickness due to north->south motion
                delHrt(2:end,:) = delH(1:end-1,:) .* dZidy_right(1:end-1,:)./gradT(1:end-1,:) ;
                % change in thickness due to south->north motion
                delHlt(1:end-1,:) = delH(2:end,:) .* dZidy_left(2:end,:)./gradT(2:end,:) ;
                
                % Track age of rockfall debris
                if LSMODE == 1
                    RFAgeup(:,2:end)   = delHup(:,2:end) .* RFAge(:,1:end-1);
                    RFAgedn(:,1:end-1) = delHdn(:,1:end-1) .* RFAge(:,2:end);
                    RFAgelt(1:end-1,:) = delHlt(1:end-1,:) .* RFAge(2:end,:);
                    RFAgert(2:end,:)   = delHrt(2:end,:) .* RFAge(1:end-1,:);
                    % age * thicnesses -> divde out below by total
                    % thickness for weighted age
                    RFAgeTmp = RFAge .* Ho + RFAgedn + RFAgeup + RFAgelt + RFAgert;
                end
                % final new thickness - again positive definite - sum of
                % original thickness with subtracted changes, plus
                % additions from north, south, east, west
                Ho = max( 0, Ho + delHdn + delHup + delHlt + delHrt ) ;
                
                if LSMODE == 1
                    topo = (topo - RFH) + Ho ;
                    RFH = Ho ;
                    RFAge(Ho > 0) = RFAgeTmp(Ho > 0) ./ Ho(Ho > 0);
                    RFAge(Ho <= 0) = t;
                else
                    topo = (topo - Regolith_H) + Ho ;                
                    Regolith_H = Ho + (Regolith_H<0).*Regolith_H ;  
                end
                
            end
            
            if avct >= avctMax; disp(['AOR iterations maxed out at ' int2str(avct) ' (Mode ' int2str(LSMODE) ').']); end
            
%             if LSMODE == 1
% 
%             else
%               
%             end