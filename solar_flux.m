% testing solar gravity as light
c =  299792458; % [m/s] speed of light
a_e2s =  0.006; % [m/s^2] acceleration of earth 2 sun
dr  = 1e3; % [m] precision
dt  = sqrt(2*dr/a_e2s);
SW = a_e2s*dr;
dt_dtp = sqrt(1 - SW / c^2);

de1 = 152.03e9; % [m] distance from sun to earth
de2 = de1 + dr; % [m] distance from sun to earth plus precision
P_sun = 3.9e26; % [W] [kg m^2/s^2 /s] energy during duration
E_sun = P_sun * dt; % [kg m^2/s^2] energy during duration
Ade1 = 4*pi*de1^2;
Ade2 = 4*pi*de2^2;
E_e1 = E_sun / (Ade1);
E_e2 = E_sun / (Ade2);
dE_e = E_e1 - E_e2;
a_e2s_est = dE_e / dr;
disp("done");
