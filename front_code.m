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

% -------------------------------------------------------------------------
% Analytical calculation (ETH calculation below):

D = [1.0e-14; 1.0e-9]; % semi-volatile first (m2/s)
eq_time = 1.3e-3;
% note - possibly need a logarithmic spacing of time steps due to quick 
% diffusion at start of simulation
no_ts = 1.0e5; 
Cstar = 1.0e-2; % effective saturation ratio of semi-volatile (ug/m3 (air))
Nm = 2.0e9; % number concentration of particles (#/m3)
Rpm = 1.1e-8; % particle radius (m)
% starting mass of aerosol (ug/m3 (air)), second from last term is 
% conversion from m3 to mass (e.g. *1e12 if density(p)=1 g/cm3) last term
% is number concentration of particles (#/m3)
ml0 = 4.0/3.0*pi*Rpm^3*1.0e12*Nm; 
Dg = 5.0e6; % gas phase diffusion coefficient (m2/s)
[Ca0store, times, Cg0store] = analyt_eq31_simple_Gatzsche(D, eq_time, no_ts, Cstar, ml0, Rpm, Dg);
hold on
plot(times, Ca0store(1, :), 'g','LineWidth',2)
xlabel('Time (s)', 'interpreter', 'latex', 'fontsize', 25)
ylabel('$\rm [Org]\quad (\mu g\, m^{-3}\quad (air))$', 'interpreter', 'latex', 'fontsize', 25)
set(gca, 'FontSize', 25, 'FontName', 'Cambria') % set font size

% ETH calculation:

ts = (eq_time)/1000; % initial time step for solution (s)
SN_part = 20;
Dp = Rpm*2.0; % particle diameter (m)
Dmethod = 2;
% use self-diffusion coefficients from above (row vector) (m^2 s^{-1})
Db = D.'; 
ut = 2.0e0;
estime = [0.0 eq_time]; % saturation ratio times (s)
% saturation ratios of semi-volatile component initially in bulk (1st
% element) and at surface (2nd element)
es0 = [0.00 Cstar/ml0];

% component molar masses (g mol^{-1}), 1st two are water and sucrose  
% from CRC online handbook
M = [1000.0 18.015];
p = [1.0e6 1.0e6]; % component densities (g/m3) (keep the same as in the analytical approach)

% ideality marker (1 for ideal, 0 for non-ideal)
idma = 1;

[time, time2, ch2, SN_part, Vrec, nrec, Csvrec] = ETH_Gatzsche3(ts, SN_part, Dp, Db, ut, ...
        es0, M, p, Nm, Cstar, 2.0);

%[Zrec, time, time2, ch2, SN_part, rcsave] = ETH_Gatzsche(ts, SN_part, Dp, Dmethod, Db, ut, ...
%        es0, idma, M, p, Nm);
hold on 
plot(time2, Csvrec, '-r','LineWidth',2)
legend('Zaveri','ETH') % include legend
axis([-0.05e-3 1.5e-3 -0.05 2.05])
legend boxoff % hides legend's box