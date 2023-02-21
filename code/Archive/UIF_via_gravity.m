function UIF_via_gravity()
% Code designed to demonstrate detecting a universal inertial reference
% frame (UIF) is possible via experiment
clear all
clc
close all
% constants
c  =  299792458; % [m/s] speed of light
G  = 6.6744e-11; % [m^3/(kg s)] gravitational constant
Me = 5.97219e24; % [kg] earth's mass
Ms = 333000*Me;  % [kg] sun's mass

% functions
gamma       = @(v) 1./sqrt(1-v.^2);
UIF_vel     = @(v1_in,v2_in) (v1_in+v2_in)/(1 + v1_in*v2_in);
rel_vel     = @(vUIF_in,v2_in) (v2_in-vUIF_in)/(v2_in*vUIF_in-1);
grav_2_dt   = @(g,r) sqrt(1-g*r/(.5*c^2));
r_2_gravSun = @(r) G*Ms/r^2;
gravimeter  = @(dtnear_dtfar,dr) (c^2/(2*dr))*(1-(dtnear_dtfar)^2);

%% experiement case 2: travel into center of sum

% set conditions
rSunE         = 152.03e9;                      % [m] distance from sun to earth
r_measure     = rSunE / 4;                     % 0.25 [AU] but in [m]
v_sun         = 0.2;                           % speed [frac of c] of sun (as seen by UIF)
v_probe1      = 0;                             % speed [frac of c] of 1st probe (as seen by UIF)
v_probe2      = UIF_vel(v_sun,v_sun-v_probe1); % speed [frac of c] of 2nd probe (as seen by UIF)
gravimeter_dr = 1000;                          % [m] how far apart clocks are when stationary (as seen by UIF)

% determine effects on gravimeter from orbit of sun
dr_sun        = gravimeter_dr/gamma(v_sun);          % how far appar (as seen by UIF)
r_f_orbit     = r_measure+dr_sun/2;                  % radial distance of clock farthest from source
r_n_orbit     = r_measure-dr_sun/2;                  % radial distance of clock nearest to source
dtn_dtf_orbit = frames_dtn_dtf(r_f_orbit,r_n_orbit); % relative time dilation between closest and farthest clock
g_m_orbit     = gravimeter(dtn_dtf_orbit,dr_sun);    % measured g from orbit

% determine effects on gravimeter from probe 1
dr_probe1        = gravimeter_dr/gamma(v_probe1);               % how far appar (as seen by UIF)
r_f_probe1       = r_measure+dr_probe1/2;                       % radial distance of clock farthest from source
r_n_probe1       = r_measure-dr_probe1/2;                       % radial distance of clock nearest to source
dtn_dtf_probe1   = frames_dtn_dtf(r_f_probe1,r_n_probe1);       % relative time dilation between closest and farthest clock
dr_correction_p1 = dr_probe1*gamma(rel_vel(v_sun,v_probe1));    % making what is thought to be a correction to dr
g_m_probe1       = gravimeter(dtn_dtf_probe1,dr_correction_p1); % measured g from probe 1

% determine effects on gravimeter from probe 2
dr_probe2        = gravimeter_dr/gamma(v_probe2);              % how far appar (as seen by UIF)
r_f_probe2       = r_measure+dr_probe2/2;                      % radial distance of clock farthest from source
r_n_probe2       = r_measure-dr_probe2/2;                      % radial distance of clock nearest to source
dtn_dtf_probe2   = frames_dtn_dtf(r_f_probe2,r_n_probe2);      % relative time dilation between closest and farthest clock
dr_correction_p2 = dr_probe2*gamma(rel_vel(v_sun,v_probe2));   % making what is thought to be a correction to dr
g_m_probe2       = gravimeter(dtn_dtf_probe2,dr_correction_p2);% measured g from probe 2

% displace results
temp = [[g_m_probe2 g_m_orbit g_m_probe1]/g_m_orbit ;...
  [1/gamma(v_probe2) 1/gamma(v_sun) 1/gamma(v_probe1)]*gamma(v_sun)];
fprintf('g(r) for probe 1: %0.2e [N/kg]\n',g_m_probe1);
fprintf('g(r) for orbit: %0.2e [N/kg]\n',g_m_orbit);
fprintf('g(r) for probe 2: %0.2e [N/kg]\n',g_m_probe2);

% supporting function
function dtn_dtf = frames_dtn_dtf(r_f,r_n)
g_f       = r_2_gravSun(r_f);        % gravitational specific force at clock farthest from source (as seen by UIF and orbit)
g_n       = r_2_gravSun(r_n);        % gravitational specific force at clock nearest to source (as seen by UIF and orbit)
dt_f      = grav_2_dt(g_f,r_f);      % time dilation of clock farthest from source (as seen by UIF)
dt_n      = grav_2_dt(g_n,r_n);      % time dilation of clock nearest to source (as seen by UIF)
dtn_dtf   = dt_n/dt_f;               % relative time dilation between closest and farthest clock
end
end
