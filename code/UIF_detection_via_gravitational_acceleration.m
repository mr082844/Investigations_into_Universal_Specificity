% Code designed to demonstrate detection of universal inertial frame (UIF)
function UIF_via_gravity()
%% initializations, constants and simple functions
% initialization
clear all
clc
close all

% constants
c  = 299792458;  % [m/s] speed of light
G  = 6.6744e-11; % [m^3/(kg s)] gravitational constant
Me = 5.97219e24; % [kg] earth's mass
Ms = 333000*Me;  % [kg] sun's mass
et = 0.5*c^2;    % [m^2/s^2] specific total energy
AU = 152.03e9;   % [m] distance from sun to earth

% simple functions
gamma       = @(v) 1./sqrt(1-v.^2);
add_vel     = @(v1_in,v2_in) (v1_in+v2_in)/(1 + v1_in*v2_in);
grav_2_dt   = @(g,r) sqrt(1-g*r/et);
r_2_gravObj = @(M,r) G*M/r^2;
gravimeter  = @(dtnear_dtfar,dr) (c^2/(2*dr))*(1-(dtnear_dtfar)^2);

%% experiement: travel two gravimeters (probes) towards center of massed object (MO)
% set conditions (in MO's frame)
M_MO      = 1e3*Ms; % [kg] mass of object at center of experiment
r_measure = AU/2;   % [m] nearest clock distance from center of MO
probe_dv  = 0.1;    % [frac of c] speed of probes relative to MO
gmtr_dr   = 1000;   % [m] clocks distance apart when stationary

% initialize
gr_orbit_all  = [];
gr_probe1_all = [];
gr_probe2_all = [];

% loop through range of MO velocities
v_obj_all = [0:0.01:0.99 0.99:0.001:0.999]; % [frac of c] speed of MO (in UIF)
for ivo = 1 : length(v_obj_all)
    % (in UIF)
    v_obj         = v_obj_all(ivo);           % [frac of c] velocity of MO
    v_p1          = add_vel(v_obj,probe_dv);  % [frac of c] velocity of probe1
    v_p2          = add_vel(v_obj,-probe_dv); % [frac of c] velocity of probe2
    drUIF_drp_obj = gamma(v_obj);             % [-] kinetic differential for MO
    drUIF_drp_p1  = gamma(v_p1);              % [-] kinetic differential for probe1
    drUIF_drp_p2  = gamma(v_p2);              % [-] kinetic differential for probe2
    
    % determine kinetic time/space dilation effects on grivimeters (in UIF)
    gmtr_dr_UIF_obj = gmtr_dr/drUIF_drp_obj; % [m] clocks distance apart
    gmtr_dr_UIF_p1  = gmtr_dr/drUIF_drp_p1;  % [m] clocks distance apart
    gmtr_dr_UIF_p2  = gmtr_dr/drUIF_drp_p2;  % [m] clocks distance apart
    
    % determine effects on gravimeter from orbit of MO (in  MO frame)
    dr_obit       = drUIF_drp_obj*gmtr_dr_UIF_obj;       % [m] clocks distance apart
    r_f_orbit     = r_measure+dr_obit;                   % [m] farthest clock distance to MO
    r_n_orbit     = r_measure;                           % [m] nearest clock distance to MO
    dtn_dtf_orbit = frames_dtn_dtf(r_f_orbit,r_n_orbit); % [-] clock differential
    g_m_orbit     = gravimeter(dtn_dtf_orbit,gmtr_dr);   % [m/s^2] measured g
    
    % determine effects on gravimeter from probe 1 (in  MO frame)
    dr_p1          = drUIF_drp_obj*gmtr_dr_UIF_p1;          % [m] clocks distance apart
    r_f_probe1     = r_measure+dr_p1;                       % [m] farthest clock distance to MO
    r_n_probe1     = r_measure;                             % [m] nearest clock distance to MO
    dtn_dtf_probe1 = frames_dtn_dtf(r_f_probe1,r_n_probe1); % [-] clock differential
    g_m_probe1     = gravimeter(dtn_dtf_probe1,gmtr_dr);    % [m/s^2] measured g
    
    % determine effects on gravimeter from probe 2 (in  MO frame)
    dr_probe2      = drUIF_drp_obj*gmtr_dr_UIF_p2;          % [m] clocks distance apart
    r_f_probe2     = r_measure+dr_probe2;                   % [m] farthest clock distance to MO
    r_n_probe2     = r_measure;                             % [m] nearest clock distance to MO
    dtn_dtf_probe2 = frames_dtn_dtf(r_f_probe2,r_n_probe2); % [-] clock differential
    g_m_probe2     = gravimeter(dtn_dtf_probe2,gmtr_dr);    % [m/s^2] measured g
    
    % store results
    gr_orbit_all  = [gr_orbit_all g_m_orbit];
    gr_probe1_all = [gr_probe1_all g_m_probe1];
    gr_probe2_all = [gr_probe2_all g_m_probe2];
end

% plot results
fig = figure(1);
hold off
plot(v_obj_all, gr_orbit_all,'LineWidth',2);
hold on
plot(v_obj_all, gr_probe1_all,'LineWidth',2);
plot(v_obj_all, gr_probe2_all,'LineWidth',2);

% clean up plot
legend('Orbital g(r)','Gravimeter1 g(r)','Gravimeter2 g(r)','FontSize',16,'location','NW');
xlabel({'Dimensional Velocity of Massed Object [fraction of c]','With Respect to UIF'},'FontSize',16);
ylabel('Measured $g(r)~\left[\frac{m}{s^2}\right]$','FontSize',16,'Interpreter','latex');
grid on
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'fontsize',16)
xticks([0:.1:1]);
annotation(fig, 'textbox', [.13 .10 .8 .2], 'String'...
    ,sprintf('Mass of massed object: %d [Solar Masses]',M_MO/Ms)...
    ,'EdgeColor','none','FontSize',14);
annotation(fig, 'textbox', [.13 .07 .8 .2], 'String'...
    ,sprintf('Gravimeter clocks original distance apart: %d [km]',gmtr_dr/1e3)...
    ,'EdgeColor','none','FontSize',14);
annotation(fig, 'textbox', [.13 .04 .8 .2], 'String'...
    ,sprintf('Measurement distance from center of mass: %0.1f [AU]',r_measure/AU)...
    ,'EdgeColor','none','FontSize',14);
annotation(fig, 'textbox', [.13 .01 .8 .2], 'String'...
    ,sprintf('Speed of probes relative to massed object: %0.1f [fraction of c]',probe_dv)...
    ,'EdgeColor','none','FontSize',14);

%% supporting function
    function dtn_dtf = frames_dtn_dtf(r_f,r_n)
        % (in MO frame)
        g_f     = r_2_gravObj(M_MO,r_f); % gravitational specific force at clock farthest from MO
        g_n     = r_2_gravObj(M_MO,r_n); % gravitational specific force at clock nearest to MO
        dt_f    = grav_2_dt(g_f,r_f);    % time dilation of clock farthest from MO
        dt_n    = grav_2_dt(g_n,r_n);    % time dilation of clock nearest to MO
        dtn_dtf = dt_n/dt_f;             % relative time differential between closest and farthest clock
    end
end