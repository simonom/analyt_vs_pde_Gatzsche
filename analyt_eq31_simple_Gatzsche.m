% Copyright Notice: This code is in Copyright.  Any use leading to 
% publication or 
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
% -------------------------------------------------------------------------

% the equivalent of the analytical solution that is nested in 
% Zaveri_Eq26Eq31_pde_D_dependent, so that a fair comparison can be made
% to the analytical solution there, with its effective D derived from the
% pde solution, rather than use analyt_eq31_D_dependent, which is a 
% slightly different code for the analytical solution because it considers
% volume change during the Newton step method for iteration

% different to analyt_eq31_simple because it does not calculate D as a
% function of composition, but it does take inputs

% this version for use with analyt_vs_pde_Zaveri_compare.m
function [Ca0store, times, Cg0store, Dstore] = ...
    analyt_eq31_simple_Gatzsche(Db, eq_time, no_ts, Cstar, ml0, Rpm, Dg)

    % inputs:
    % D - self-diffusion coefficients of components (semi-volatile first)
    % eq_time - diffusion equilibrium time (s)
    % no_ts - number of time steps to solve diffusion over
    % Cstar - effective saturation concentration of semi-volatile (ug/m^3 
    % (air))
    % ml0 - initial mass loading of condensed phase (ug/m^3 (air))
    
    % starting values and constants
    % result times (s)
    times = horzcat(0,logspace(log10(eq_time)*4.0, log10(eq_time), no_ts)); 
    %times = linspace(0, eq_time, no_ts); 
    % starting values and constants
    [kci, Db, Rpm, p, Cstar, Ai, Cg, Nm, ~, M] = ...
        starting_stuff(times, Db, Cstar, ml0, Rpm);
    % x mesh (m) and shell volumes (m^3)
    [Vpt] = spat_arrays(Rpm, Nm);
    
    
    % records of important values used by the analytical solution through
    % time steps.  Dbstore has a length one less than the times because the
    % Db values represent those used over each time interval, rather than
    % those at the start and end of each time interval.
    Vp0store = zeros(1, length(times));
    Rpm0store = zeros(1, length(times));
    Cg0store = zeros(1, length(times));
    Dstore = zeros(1, length(times));
    % particle-phase concentration (ug/m^3 (air))
    Ca0store = zeros(2, length(times)); 
    % remember initial values
    Vp0store(1) = Vpt(1);
    Rpm0store(1) = Rpm;
    Cg0store(1) = Cg(1, 1);
    % particle-phase concentration at t=0 s (ug/m^3 (air))
    Ca0store(:, 1) = Ai(:, 1).*(Vpt(1));
    
    % sum of the concentration of all components in the particle-phase 
    % (ug/m^3 (particle))
    Aj = p; % when all components have the same density
    % mass of non-volatile component in a single particle 
    % (ug/(single particle))
    mugnv = Ca0store(2, 1)/Nm;
    
    for ti = 1:length(times)-1
        
        disp(ti)
        h = times(ti+1)-times(ti); % time step (s)
        % recall the total particle-phase volume and the radius from the 
        % end of the previous time step (m^3)
        Vp0 = Vp0store(ti); 
        Rpm0 = Rpm0store(ti);
        Cg0 = Cg0store(ti); Cg2 = Cg0store(ti);
        Ca0 = Ca0store(:, ti); Ca2 = Ca0store(:, ti);
        D = Db(1)^(Ca0(1)/sum(Ca0(:)))*Db(2)^(Ca0(2)/sum(Ca0(:)));
        Dstore(ti) = D;
        % overall gas-side mass transfer coefficient (m/s) (eq. 20)
        Kg = Kgcalc(Dg, Rpm0, D, Cstar, p, kci);
        % concentration difference factor in eqs. 29/31 (/s)
        k2 = Vp0*(3.0/Rpm0)*Kg;
        
        % total number of moles in bulk particle-phase (mol)
        xT = (Ca2(1)/M(1)+Ca2(2)/M(2));
        % semi-volatile mole fraction in bulk condensed-phase (fraction)
        xsv = (Ca2(1)/M(1))/xT;
        
        % eqs. 31/32 for mass-transfer
        % average concentration in the particle-phase (ug/m} (air))
        errora = k2*(Cg2-xsv*Cstar)-(kci*Ca2(1))-...
            ((Ca2(1)-Ca0(1))/h);
        
        errorad = k2*(-(Cstar/(Aj*Vp0)))-kci-1/h; % differentiated form
        
        errorg = -k2*(Cg2-xsv*Cstar)-...
            ((Cg2-Cg0)/h);%+gamma*(h/3.6e3);
        errorgd = -k2-1/h; % differentiated form
        
        % loop until we get the error acceptably low - note error tolerance
        % may need to vary depending on relative diffusion coefficients of
        % components
        while abs(errora)>(Ca2(1)/1.0e0) || abs(errorg)>(Cg2/1.0e0)

            % improve concentration at end of time step estimates (ug/m3 
            % (air))
            Ca2(1) = Ca2(1)-errora/errorad;
            Cg2(1) = Cg2(1)-errorg/errorgd;

            % total number of moles in condensed-phase (mol)
            xT = (Ca2(1)/M(1)+Ca2(2)/M(2));
            % semi-volatile mole fraction in condensed-phase (fraction)
            xsv = (Ca2(1)/M(1))/xT;
            
            % revalue error
            errora = k2*(Cg2-xsv*Cstar)-(kci*Ca2(1))-...
            ((Ca2(1)-Ca0(1))/h);
        
            errorg = -k2*(Cg2-xsv*Cstar)-...
            ((Cg2-Cg0)/h);%+gamma*(h/3.6e3);
            errorgd = -k2-1/h; % differentiated form
     
        end
        
        % new mass of semi-volatile (ug/single particle)
        mugsv = Ca2(1)/Nm;
        
        % new volume of single particle (m^3)
        Vsingle = (mugsv+mugnv)*(1/p);
%         Vsingle = Vp0/Nm;
        % new particle bulk-average concentrations (ug/m^3 (particle))
        C1 = mugsv/Vsingle;
        Cnv = mugnv/Vsingle;
        
        % new volume of all particles (m^3 m^{-3} (air))
        Vp0store(ti+1) = Vsingle*Nm;
        % new concentration of condensed sv and nv (ug/m^3 (air))
        Ca2(1) = C1*Vp0store(ti+1);
        Ca2(2) = Cnv*Vp0store(ti+1);

        % new particle radius (m)
        Rpm0store(ti+1) = ((3*Vp0store(ti+1))/(4*pi*Nm))^(1/3);
        % store particle-phase concentration (ug/m^3 (air))
        Ca0store(:, ti+1) = Ca2;
        % store gas-phase concentration (ug/m^3 (air)
%         Cg0store(ti+1) = Cg0store(ti); % infinite mass of sv
        % finite mass of semi-volatile (ensure gas-phase concentration
        % estimation is turned on in the iteration part)
        Cg0store(ti+1) = Cg0store(ti)-(Ca0store(1, ti+1)-Ca0store(1,ti)); 
        
        % finite mass of sv with emission (ug/m^3 (air)
%         Cg0store(ti+1) = Cg0store(ti)-...
%             (Ca0store(1,ti+1)-Ca0store(1,ti))+gamma*(h/3.6e3);
        
        
    end % time loop

end % analytical function

    
% -------------------------------------------------------------------------
% nested function containing constants for input to the partitioning
% equation
function [kci, D0, Rpm, p, Cstar, Ai, Cg, Nm, gamma, M] = ...
    starting_stuff(times, D, Cstar, ml0, Rpm)

    % inputs:
    % times - time points to solve diffusion at (s)
    % D - self-diffusion coefficients of components,semi-volatile first
    % (ug/m^3 (air))
    % Cstar - effective saturation concentration of semi-volatile (ug/m^3 
    % (air))
    
    % particle-phase bulk concentration (ug/m^3 (particle)) and gas-phase 
    % bulk concentrations (ug/m^3 (air)) 
    Ai = zeros(2, length(times));
    Cg = zeros(2, length(times));

    p = 1.0e12; % component densities (assumed equal) (ug m^{-3})
    M = [1.0e9 18.015e6]; % molecular weight (ug mol^{-1}) (sv first)
    
    % concentration of particles (#/m^3)
    Nm = round((3*ml0)/(4*pi*p*(Rpm^3)));
    
    kci = 0.0; % pseudo-first-order reaction constant (/s)
    % self-diffusion coefficients of components, semi-volatile first, then
    % non-volatile (m^2/s)
    D0 = D; 

    % initial gas- and particle-phase bulk concentrations (ug m^{-3}) for the
    % semi-volatile in the first row and non-volatile in the second row
    Cg(:,1) = [2.0; 0.0]; 
    Ai(:,1) = [0.0; 1.0e12]; 
    
    gamma = 0.0; % emission rate of semi-volatile (ug/m3(air)/hr)
    
end
% -------------------------------------------------------------------------
% nested function for total particle volume
function [Vpt] = spat_arrays(Rpm,Nm)
    
    % total particle-phase volume (m^3/m^3(air))
    Vpt(1) = ((4.0/3.0)*pi*Rpm^3.0*Nm); 
end
% -------------------------------------------------------------------------
% nested function to calculate gas-side mass transfer coefficient (eq. 20 
% of Zaveri et al. 2014) 
function Kg = Kgcalc(Dg, Rp, Db, Cstar, Aj, kci)

    if kci == 0     
        % particle-side mass transfer coefficient (m/s) (eq. 24)
        kp = 5.0*(Db/Rp);     
    else
        qi = Rp*((kci/Db)^0.5);
        [Qi] = Qcalc(kci, Db, Rp);
        % particle-side mass transfer coefficient (m/s)  (eq. 23)
        kp = (Db/Rp)*((qi*coth(qi)-1.0)/(1.0-Qi)); 
    end
    % *********************************************************************
    % gas-side mass transfer coefficient (m/s)
    
    % accommodation coefficient (pp. 5 of Zaveri et al. (2008)) (-)
    alpha = 1.0;
    % gas-phase mean free path (what MOSAIC gets when iv==12 in line 9147 
    % of mosaic_box.25.f90) (m)
    lambda = 5.972e-8;  
    Kn = lambda/Rp; % Knudsen number (-)
    % transition regime correction factor
    f = (0.75*alpha*(1+Kn))/(Kn*(1+Kn)+0.283*alpha*Kn+0.75*alpha);

    kg = (Dg/Rp)*f; % gas-side mass transfer coefficient (m/s)
    
    % ********************************************************************* 
    SR = (Cstar/(Aj)); % coefficient for kp term (eq. 20)
    % overall gas-side mass transfer coefficient (m/s) (eq. 20)
    Kg = 1/(1/(kg)+(1/(kp)*SR)); 
    %Kg = 1/(0+(1/(kp)*SR)); 
end