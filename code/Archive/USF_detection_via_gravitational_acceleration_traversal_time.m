% Code designed to demonstrate detection of universally stationary frame (USF)
function USF_via_gravity()
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
M_MO        = 1e3*Ms; % [kg] mass of object at center of experiment
r_measure   = AU/2;     % [m] nearest clock distance from center of MO
probe_dv    = 0.1;   % [frac of c] speed of probes relative to MO
dr_orb_MO   = 1e3;    % [m] clocks distance apart when stationary

% initialize
gr_orbit_all  = [];
gr_probe1_all = [];
gr_probe2_all = [];

% loop through range of MO velocities
v_obj_all = [0:0.01:0.99 0.99:0.001:0.999]; % [frac of c] speed of MO (in USF)
for ivo = 1 : length(v_obj_all)
    % (in USF)
    v_obj          = v_obj_all(ivo);            % [frac of c] velocity of MO
    v_p1           = add_vel(v_obj,probe_dv);   % [frac of c] velocity of probe1
    v_p2           = add_vel(v_obj,-probe_dv);  % [frac of c] velocity of probe2
    drUSF_drp_obj  = gamma(v_obj);              % [-] kinetic differential for MO
    drUSF_drp_p1   = gamma(v_p1);               % [-] kinetic differential for probe1
    drUSF_drp_p2   = gamma(v_p2);               % [-] kinetic differential for probe2
    drp_p1_drp_obj = drUSF_drp_obj/drUSF_drp_p1;% [-] kinetic differential WRT MO
    drp_p2_drp_obj = drUSF_drp_obj/drUSF_drp_p2;% [-] kinetic differential WRT MO
    
    % determine miscalibration effects on grivimeters (in MO frame)
    dr_p1_MO  = dr_orb_MO*drp_p1_drp_obj; % [m] clocks distance apart
    dr_p2_MO  = dr_orb_MO*drp_p2_drp_obj; % [m] clocks distance apart
    
    % determine distance traversal effects (in MO frame)
    c_fn_p1    = c*drp_p1_drp_obj^-2*(1-v_p1);     % miscalibrated effective speed of light
    c_fn_p2    = c*drp_p1_drp_obj^-2*(1+v_p2);     % miscalibrated effective speed of light
    dt_add_p1  = dr_p1_MO/c_fn_p1;                 % traversal time (small compared to integration)
    dt_add_p2  = dr_p2_MO/c_fn_p2;                 % traversal time (small compared to integration)
    dr_add_orb = 0;
    dr_add_p1  = dt_add_p1*probe_dv*c;     % traversal additional distance (constant)
    dr_add_p2  = dt_add_p2*probe_dv*c;     % traversal additional distance (constant)
    
    % determine effects on gravimeter from orbit of MO (in  MO frame)
    r_f_orbit     = r_measure+dr_orb_MO+dr_add_orb;      % [m] farthest clock distance to MO
    r_n_orbit     = r_measure;                           % [m] nearest clock distance to MO
    dtn_dtf_orbit = frames_dtn_dtf(r_f_orbit,r_n_orbit); % [-] clock differential
    g_m_orbit     = gravimeter(dtn_dtf_orbit,dr_orb_MO); % [m/s^2] measured g
    
    % determine effects on gravimeter from probe 1 (in  MO frame)
    r_f_probe1     = r_measure+dr_p1_MO+dr_add_p1;          % [m] farthest clock distance to MO
    r_n_probe1     = r_measure;                             % [m] nearest clock distance to MO
    dtn_dtf_probe1 = frames_dtn_dtf(r_f_probe1,r_n_probe1); % [-] clock differential
    g_m_probe1     = gravimeter(dtn_dtf_probe1,dr_orb_MO);  % [m/s^2] measured g
    
    % determine effects on gravimeter from probe 2 (in  MO frame)
    r_f_probe2     = r_measure+dr_p2_MO+dr_add_p2;          % [m] farthest clock distance to MO
    r_n_probe2     = r_measure;                             % [m] nearest clock distance to MO
    dtn_dtf_probe2 = frames_dtn_dtf(r_f_probe2,r_n_probe2); % [-] clock differential
    g_m_probe2     = gravimeter(dtn_dtf_probe2,dr_orb_MO);  % [m/s^2] measured g
    
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
xlabel({'Dimensional Velocity of Massed Object [fraction of c]','With Respect to USF'},'FontSize',16);
ylabel('Measured $g(r)~\left[\frac{m}{s^2}\right]$','FontSize',16,'Interpreter','latex');
grid on
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'fontsize',16)
xticks([0:.1:1]);
annotation(fig, 'textbox', [.13 .10 .8 .2], 'String'...
    ,sprintf('Mass of massed object: %d [Solar Masses]',M_MO/Ms)...
    ,'EdgeColor','none','FontSize',14);
annotation(fig, 'textbox', [.13 .07 .8 .2], 'String'...
    ,sprintf('Gravimeter clocks original distance apart: %d [km]',dr_orb_MO/1e3)...
    ,'EdgeColor','none','FontSize',14);
annotation(fig, 'textbox', [.13 .04 .8 .2], 'String'...
    ,sprintf('Measurement distance from center of mass: %0.1f [AU]',r_measure/AU)...
    ,'EdgeColor','none','FontSize',14);
annotation(fig, 'textbox', [.13 .01 .8 .2], 'String'...
    ,sprintf('Speed of probes relative to massed object: %0.1f [fraction of c]',probe_dv)...
    ,'EdgeColor','none','FontSize',14);
set(gca,'YScale','log')
ylim([10^0 10^4]);

%% supporting function
    function dtn_dtf = frames_dtn_dtf(r_f,r_n)
        % (in MO frame)
        g_f       = r_2_gravObj(M_MO,r_f); % gravitational specific force at clock farthest from MO
        g_n       = r_2_gravObj(M_MO,r_n); % gravitational specific force at clock nearest to MO
        dt_f      = grav_2_dt(g_f,r_f);    % time dilation of clock farthest from MO
        dt_n      = grav_2_dt(g_n,r_n);    % time dilation of clock nearest to MO
        dtn_dtf   = dt_n/(dt_f);    % relative time differential between closest and farthest clock
    end
end