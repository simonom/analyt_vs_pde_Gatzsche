% Copyright Notice: This code is in Copyright.  Any use leading to publication or 
% financial gain is prohibited without the permission of the authors Simon 
% O'Meara : simon.omeara@manchester.ac.uk.  First published 2017.

% This file is part of analyt_vs_pde_Gatzsche

% analyt_vs_pde_Gatzsche 
% is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.

% analyt_vs_pde_Gatzsche
% is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.

% You should have received a copy of the GNU General Public License
% along with analyt_vs_pde_Gatzsche 
% (see the LICENSE file).  If not, see 
% <http://www.gnu.org/licenses/>.

% Code for checking Gatzsche assumption that Zaveri equation with 
% component-dependent D 
% equal to numerical solution with composition dependent D 
% Double precision on inputs to the diffusion equation.

% ******************DOUBLE PRECISION VERSION*******************************

% Simon O'Meara, Unviersity of Manchester, May 2015.

function [time, time2, check2, SN_part, Vrec, nrec, Csvrec] = ETH_Gatzsche3(ts, SN_part, Dp, Db, ut, es0, M, p, Nm, Cstar, Cg0)

    % ---------------------------------------------------------------------
    % constants inputted at the command line:

    % ts - time step (s)
    % SN_part - number of shells
    % Dp - initial particle diameter (cm)
    % Dmethod - diffusion coefficient dependence on mole fraction
    % Db - self-diffusion coefficients of the semi-volatile (first) and  
    % non-volatile (latter) (m^2 s^{-1}) 
    % ut - acceptable percentage change in number of moles in a single
    % shell per time step
    % RH - relative humidity (fraction) in first row
    % idealmark - determines whether to use equations for ideality (1) or
    % non-ideality (0)
    
    time = zeros(1, 1); % time array (s) (grows with time steps)
     
    sp = 1.0e3; % number of time points we want to save results at
    % time array to relate to Zrec array (s).
    time2 = zeros(1, sp);
    
    % check for equation breaking
    check2 = 0;
    % display option
    format long
    
    % ---------------------------------------------------------------------
    % prepare arrays
    
    % spatial arrays
    [delta0, ~, A, V0, ~] = spat_arrays(SN_part, Dp, Cg0, Cstar);
    % initial mole fraction per shell (columns) and component (rows (sv 1st))
    xi0 = ones(2, SN_part);
    xi0(1, :) = horzcat(xi0(1, 1:SN_part-1)*es0(1), Cg0/Cstar);
    
    if Cg0/Cstar >= 0.0 % supersaturated case
        xi0(2, :) = horzcat(xi0(2, 1:SN_part-1)*(1.0-es0(1)), 0.0);
    end
    
    % initial concentration and number of moles array (mol m^{-3) and mol),
    % also initial diffusion coefficient array
    [Z0, Nz0, nrec, Vrec, Csvrec] = concs(xi0, M, SN_part, V0, sp, p, Nm, Cg0);
    
    % ---------------------------------------------------------------------
    % time loop
    tsn = 0; % count on time steps
    ts0 = ts; % original time step interval (s)
    zreci = 1; % count on result storing
    %----------------------------------------------------------------------
    % time loop
    while time(tsn+1) < (ts0*sp)
        
        % update time array (s)
        tsn = tsn+1;
        time(tsn+1) = time(tsn)+ts;
        
        % new saturation ratio of semi-volatile
        Z1 = Z0; Nz1 = Nz0; delta1 = delta0;
        % mole fraction in each shell
        xsv = Nz1(1, :)./sum(Nz1, 1);
        xnv = Nz1(2, :)./sum(Nz1, 1);
        % mole fraction at shell boundaries
        xsv = (xsv(2:SN_part)+xsv(1:SN_part-1))/2.0;
        xnv = (xnv(2:SN_part)+xnv(1:SN_part-1))/2.0;
        % diffusion coefficient at shell boundaries (m^2 s^{-1})
        Db1 = Db(1).^xsv.*Db(2).^xnv;
        
        % outputs absolute number of moles and total volume per shell
        [Nznew, Vnew, Znew] = eqn(A, SN_part, M, Nz1, delta1, Db1, Z1, ts, p);
        
        %------------------------------------------------------------------
        % percentage change in number of moles of sv per shell 
        % over the current time step
        x_chsv = abs(Nznew(1, :)./sum(Nznew(:, :), 1)-Nz1(1, :)./sum(Nz1(:, :), 1));
        x_chnv = abs(Nznew(2, :)./sum(Nznew(:, :), 1)-Nz1(2, :)./sum(Nz1(:, :), 1));
        % find where change is exceeded
        ex_isv = x_chsv>ut;
        ex_inv = x_chnv>ut;
        
        % if acceptable change is not exceeded and we're not at the maximum
        % time step then try increasing time step.       
        while sum(ex_isv) == 0 && sum(ex_inv) == 0 && ts<ts0
            
            ts = ts*2.0;
            time(tsn+1) = time(tsn)+ts;

            % outputs absolute number of moles and total volume per shell
            [Nznew, Vnew, Znew] = eqn(A, SN_part, M, Nz1, delta1, Db1, Z1, ts, p);

            % percentage change in number of moles of sv per shell 
            % over the current time step
            x_chsv = abs(Nznew(1, :)./sum(Nznew(:, :), 1)-Nz1(1, :)./sum(Nz1(:, :), 1));
            x_chnv = abs(Nznew(2, :)./sum(Nznew(:, :), 1)-Nz1(2, :)./sum(Nz1(:, :), 1));

            % find where change is exceeded
            ex_isv = x_chsv>ut;
            ex_inv = x_chnv>ut;
            
        end
        
        
        % if acceptable change is exceeded decrease time step
        while sum(ex_isv)>0 || sum(ex_inv)>0 || sum(sum(Nznew<0))>0

            ts = ts/2.0;
            time(tsn+1) = time(tsn)+ts; 
            
            % outputs absolute number of moles and total volume per shell
            [Nznew, Vnew, Znew] = eqn(A, SN_part, M, Nz1, delta1, Db1, Z1, ts, p);

            % percentage change in number of moles of sv per shell 
            % over the current time step
            x_chsv = abs(Nznew(1, :)./sum(Nznew(:, :), 1)-Nz1(1, :)./sum(Nz1(:, :), 1));
            x_chnv = abs(Nznew(2, :)./sum(Nznew(:, :), 1)-Nz1(2, :)./sum(Nz1(:, :), 1));
            % find where change is exceeded
            ex_isv = x_chsv>ut;
            ex_inv = x_chnv>ut;
            
        end
        
        
        Nz0 = Nznew; % revalue number of moles (mol)
        V0 = Vnew; % revalue shell volumes (m^3)
        % revalue concentrations at start of time step (mol m^{-3})
        Z0 = Znew;
        
        % -----------------------------------------------------------------
        % revalue spatial variables
        
        % new spatial arrays.
        [delta0, ~, A] = spat_arrays_inequal(V0, SN_part);
        
        % update Cg
        Cg = Cg0-sum(Nz0(1, 1:SN_part-1))*Nm*M(1)*1.0e6;
     
        % new concentration in the surface shell:
          
        % volume of surface shell (m3) (arbitrary)
        V0(SN_part) = sum(V0(1:SN_part-1))*1.0e-2;
        % number of moles sv in surface shell
        if Cg/Cstar >= 1.0
            Nz0(1, SN_part) = (p(1)./M(1)).*V0(SN_part).*(Cg/Cstar);
        else
            xsv = Cg/Cstar;
            Nz0(1, SN_part) = V0(SN_part)*(1.0/(M(1)/p(1)+(((1.0-xsv)/xsv)*M(2)/p(2))));
        end
        Z0(1, SN_part) = Nz0(1, SN_part)/V0(SN_part);
        

        % if the near-surface shell has grown too large make a new 
        % near-surface shell
        if  delta0(SN_part-1) > sum(delta0(1:SN_part-1))/SN_part;
            
            % volume of the new near-surface shell (m^3)
            Vnsl = V0(SN_part-1)*1.0e-2;
            % new shell volume array (m^3)
            V0 = horzcat(V0(1:SN_part-2), V0(SN_part-1)-Vnsl, Vnsl, V0(SN_part));
            % new concentration array (mol m^{-3})
            Z0 = horzcat(Z0(:, 1:SN_part-1), Z0(:, SN_part-1), Z0(:, SN_part));
            % increase number of shells
            SN_part = SN_part+1;
            % new number of moles (mol), organic first
            Nz0 = zeros(size(Z0, 1), SN_part);
            
            for ic1 = 1:size(Z0, 1) % component loop
                Nz0(ic1, :) = Z0(ic1, :).*V0;
            end
            clear ic1
            [delta0, ~, A] = spat_arrays_inequal(V0, SN_part);
            disp('new shell')
        end
        
        % record important information
        if time(tsn+1) > ts0*zreci
            zreci = zreci+1;
            nrec(zreci, :, 1:SN_part) = Nz0;
%             Vrec(zreci, 1:SN_part) = V0; 
            Vrec(zreci, 1:length(Db1)) = Db1;
            % mass concentration of sv in particle bulk (ug m^{-3} (air))
            Csvrec(zreci, 1) = sum(nrec(zreci, 1, 1:SN_part-1), 3)*Nm*M(1)*1.0e6;
            time2(1, zreci) = time(tsn+1);
            disp('zreci = ')
            disp(zreci)
            disp(sum(Nz0(1, 1:SN_part-1), 2)/sum(sum(Nz0(:, 1:SN_part-1), 2), 1))
            
        end
        
        
    end % time loop
    
    
    % ---------------------------------------------------------------------
    % nested function to revalue spatial arrays
    function [delta, rc, A] = spat_arrays_inequal(V0, SN_part)
        % all shells
        A = zeros(1, SN_part); % shell areas (m^2)
        rc = zeros(1, SN_part); % shell end points (m)
        delta = zeros(1, SN_part); % shell widths (m)
        for ir = 1:SN_part
            delta(1, ir) = ((3.0*sum(V0(1:ir)))/(4.0*pi))^(1.0/3.0)-...
                ((3.0*sum(V0(1:ir-1)))/(4.0*pi))^(1.0/3.0);
            rc(1, ir) = ((3.0*sum(V0(1:ir)))/(4.0*pi))^(1.0/3.0);
            A(1, ir) = 4.0*pi*rc(1, ir)^2.0;
        end
        clear ir
        delta(SN_part) = sum(delta(1:SN_part-1))/2.0; % width of particle phase film (m)
    end
    % ---------------------------------------------------------------------
    % nested function to value initial spatial arrays
    function [delta, rc, A, V0, Diamw] = spat_arrays(SN_part, Dp, Cg, Cstar)
        
        % width of surface shell (m)
        Diamw = (Dp/2.0)/1.0e3;
        r0 = (Dp/2.0)-(Diamw); % initial bulk radius (m)
        rc = zeros(SN_part, 1); % cumulative shell widths (m)
        delta = zeros(SN_part, 1); % radial width of shells (m)

        for ir2 = 1:SN_part-1
            rc(ir2) = (r0/(SN_part-1))*ir2;
            delta(ir2) = (r0/(SN_part-1));
        end
        clear ir2
        
        
        rc(end) = (Dp/2.0);
        
        if Cg/Cstar >= 0.0 % supersaturated case
            % width of particle-phase film (m)
            delta(end) = r0/2.0; 
        else
            delta(end) = Diamw;
        end
        % areas of shells (m^2)
        A = (4.0*pi*(rc.^2.0)).';
        % volume of individual shells (m^3)
        V0 = horzcat((4.0/3.0)*pi*(rc(1).')^3.0,((4.0/3.0)*pi)*((rc(2:end).').^3.0 ...
            -(rc(1:length(rc)-1).').^3.0));
    end
   
    % ---------------------------------------------------------------------
    % nested function for calculating the initial concentration
    % in the particle phase (mol m^{-3})
    function [Z, Nz0, nrec, Vrec, Csvrec] =  concs(xi, M, SN_part, V0, sp, p, Nm, Cg0)
        
        % **************************
        % inputs 
        % xi - mole fractions of each component at start 
        % with components in 1st dim. and shells in 2nd dim.
        % M - molar mass of each component (g/mol) (row array)
        % p - density of each component (g/m^3) (row array)
        % idma - ideality marker (1 for ideal, 0 for non-ideal)
        % **************************
        % outputs
        % Z - mol m^{-3} (per component (1st dim.), per shell (2nd dim.))
        % Nz0 - absolute number of moles per component (1st dim.), per
        %       shell (2nd dim.)
        % nrec - record of component concentrations (1st dim.), with shell
        %        (2nd dim.) and with time (3rd dim.) (mol m^{-3})
        % Vrec - record of shell volumes (1st dim.), with time (2nd dim.)
        %        (m^{3})
        % **************************
        
        % concentration matrix (mol m^{-3})
        Z = zeros(length(M), SN_part);
        % absolute number of moles per component, per bulk shell (mol)
        Nz0 = zeros(length(M), SN_part);
        
        % shell loop
        for ir = 1:SN_part
            
            % ratios of component molar volumes ((M/p) is m^3 mol^{-1})
            mVrat = ...
                (xi(:, ir).*((M./p).'))./(sum((xi(:, ir).*((M./p).'))));
            % 
            
            % volumes of each component (m^{3})
            Vi = mVrat.*V0(ir);
            % number of moles per component at start (n) (mol) (sv in
            % 1st row)
            Nz0(:, ir) = ((p./M).').*Vi;
            % for supersaturated case
            if ir == SN_part && xi(1, SN_part) >= 1.0
                Nz0(:, ir) = ((p./M).').*Vi.*xi(:, ir);
            elseif ir == SN_part && xi(1, SN_part) < 1.0
                Nz0(1, ir) = V0(ir)*(1.0/(M(1)/p(1)+(((1.0-xi(1, ir))/xi(1, ir))*M(2)/p(2))));
                Nz0(2, ir) = 0.0; 
            end
            % concentration of each component at start (mol m^{-3})
            Z(:, ir) = Nz0(:, ir)./V0(ir);
                
            
        end
        clear ir mVrat Vi
        
        % save initial concentrations (mol m^{-3}) (components 1st dim.,
        % number of shells 2nd dim., time steps 3rd dim.)
        nrec = ones(sp, length(M), SN_part*3).*NaN; 
        nrec(1, :, 1:SN_part) = Nz0;
        
        % save initial shell volumes (shells 1st dim., time steps 3rd 
        % dim.)
        Vrec = ones(sp, SN_part*3).*NaN; 
        Vrec(1, 1:SN_part) = V0;
        % matrix for saving mass concentration of sv (ug m^{-3} (air))
        Csvrec = zeros(sp, 1);
        Csvrec(1, 1) = sum(Nz0(1, 1:SN_part-1))*Nm*M(1)*1.0e6;
        % matrix for saving gas phase concentration of sv (ug m^{-3} (air))
        Cgstore = zeros(sp, 1);
        Cgstore(1, 1) = Cg0;
            
    end

    % ---------------------------------------------------------------------
    % Nested flux equation, for eq. 7 of Zobrist et al. (2011).  Note,
    % only water diffuses in their model.
    function [Nznew, V3, Z] = eqn(A, SN_part, M, Nz1, delta, Db, Z, ts, p)
        
        % empty arrays for new number of moles per shell, and flux
        Nzsv2 = zeros(1, SN_part);
        svflux = zeros(1, SN_part);
        Nznv2 = zeros(1, SN_part);
        nvflux = zeros(1, SN_part);
        % new volumes of shells (m^3)
        V3 = zeros(1, SN_part);
        
        % loop through shell boundaries to find the flux across them
        for ir4 = 1:SN_part % shell loop
            % for bulk shells
            if ir4<SN_part
                % flux across boundaries (mol s^{-1})
                svflux(1, ir4+1) = (A(ir4))*Db(ir4)*((Z(1, ir4+1)-Z(1, ir4))/(0.5*(delta(ir4+1)+delta(ir4))));
                nvflux(1, ir4+1) = (A(ir4))*Db(ir4)*((Z(2, ir4+1)-Z(2, ir4))/(0.5*(delta(ir4+1)+delta(ir4))));
                
                % prevent nonvolatile from moving outward into surface
                if ir4 == SN_part-1
                    nvflux(1, ir4+1) = 0.0;
                end
                
                % new number of moles in a shell (sv in top row, nv in bottom) (mol)
                Nzsv2(1, ir4) = Nz1(1, ir4)+((svflux(1, ir4+1)-svflux(1, ir4)))*ts;
                Nznv2(1, ir4) = Nz1(2, ir4)+((nvflux(1, ir4+1)-nvflux(1, ir4)))*ts;
                V3(1, ir4) = (Nzsv2(1, ir4))*(M(1)/(p(1)))+(Nznv2(1, ir4))*(M(2)/(p(2)));
                
            else % for surface shell
                % new number of moles in a shell (sv in top row, nv in bottom) (mol)
                Nzsv2(1, ir4) = Nz1(1, ir4)+((-svflux(1, ir4)))*ts;
                Nznv2(1, ir4) = Nz1(2, ir4)+((-nvflux(1, ir4)))*ts;
                V3(1, ir4) = (Nzsv2(1, ir4))*(M(1)/(p(1)))+(Nznv2(1, ir4))*(M(2)/(p(2)));
            end
        end
        clear ir4
        svflux = sum(svflux(1, 1:SN_part).*ts, 2);
        nvflux = sum(nvflux, 2);
        % concatenate the organic and water number of moles arrays
        Nznew = vertcat(Nzsv2, Nznv2);
        clear Nzsv2 Nznv2
        
        % new concentrations (mol m^{-3})
        Z = zeros(2, SN_part);
        for ic = 1:length(M)
            Z(ic, :) = Nznew(ic, :)./V3;
        end

    end
    % ---------------------------------------------------------------------
end