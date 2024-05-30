% Code designed to demonstrate detection of universally stationary frame (USF)
function USF_detection_via_gravitational_acceleration()
%% initializations, constants and simple functions
% initialization
clear all
clc
close all

% constants
c   = 299792458;  % [m/s] speed of light
G   = 6.6744e-11; % [m^3/(kg s)] gravitational constant
Me  = 5.97219e24; % [kg] earth's mass
Ms  = 333000*Me;  % [kg] sun's mass
e_T = 0.5*c^2;    % [m^2/s^2] specific total energy
AU  = 152.03e9;   % [m] distance from sun to earth

% simple functions
r_s         = @(M) G*M/e_T;
gamma_inv_K = @(v) 1./sqrt(1-v.^2);
add_vel     = @(v1_in,v2_in) (v1_in+v2_in)/(1 + v1_in*v2_in);
r_2_gravObj = @(M,r) G*M/r^2;
gamma_inv_P = @(M,r) sqrt(1-r_s(M)./r);
gravimeter  = @(dtnear_dtfar,dr) (e_T/dr)*(1-(dtnear_dtfar)^2);
prop_dist   = @(M,r1,r2) (r2*gamma_inv_P(M,r2) + 0.5*r_s(M)*log(2*r2*(gamma_inv_P(M,r2)+1)-r_s(M))) ...
    - (r1*gamma_inv_P(M,r1) + 0.5*r_s(M)*log(2*r1*(gamma_inv_P(M,r1)+1)-r_s(M)));

%% experiement: travel two gravimeters (probes) towards center of massed object (MO)
% set conditions (in MO's frame)
M_MO        = 1e3*Ms; % [kg] mass of object at center of experiment
r_measure   = AU/2;   % [m] nearest clock distance from center of MO
probe_dv    = 0.1;    % [frac of c] speed of probes relative to MO
dr_orb_MO_0 = 1e3;    % [m] clocks distance apart when stationary in zero gravity

% initialize
gr_orbit_all  = [];
gr_probe1_all = [];
gr_probe2_all = [];

% loop through range of MO velocities
v_obj_all = [0:0.01:0.99 0.99:0.001:0.999]; % [frac of c] speed of MO (in USF)
for ivo = 1 : length(v_obj_all)
    % (in USF)
    v_obj          = v_obj_all(ivo);             % [frac of c] velocity of MO
    v_p1           = add_vel(v_obj,probe_dv);    % [frac of c] velocity of probe1
    v_p2           = add_vel(v_obj,-probe_dv);   % [frac of c] velocity of probe2
    drUSF_drp_obj  = gamma_inv_K(v_obj);         % [-] kinetic differential for MO
    drUSF_drp_p1   = gamma_inv_K(v_p1);          % [-] kinetic differential for probe1
    drUSF_drp_p2   = gamma_inv_K(v_p2);          % [-] kinetic differential for probe2
    drp_p1_drp_obj = drUSF_drp_obj/drUSF_drp_p1; % [-] kinetic differential WRT MO
    drp_p2_drp_obj = drUSF_drp_obj/drUSF_drp_p2; % [-] kinetic differential WRT MO

    % determine miscalibration effects on grivimeters (in MO frame)
    r2        = solve_for_r2(M_MO,r_measure,dr_orb_MO_0); 
    dr_orb_MO = (r2-r_measure);           % [m] clocks distance apart in gravity, no velocity
    dr_p1_MO  = dr_orb_MO*drp_p1_drp_obj; % [m] clocks distance apart in gravity with velocity
    dr_p2_MO  = dr_orb_MO*drp_p2_drp_obj; % [m] clocks distance apart in gravity with velocity
    
    % determine effects on gravimeter from orbit of MO (in  MO frame)
    r_f_orbit     = r_measure+dr_orb_MO;                      % [m] farthest clock distance to MO
    r_n_orbit     = r_measure;                                % [m] nearest clock distance to MO
    dtn_dtf_orbit = frames_dtn_dtf(M_MO,r_f_orbit,r_n_orbit); % [-] clock differential
    g_m_orbit     = gravimeter(dtn_dtf_orbit,dr_orb_MO);      % [m/s^2] measured g
    
    % determine effects on gravimeter from probe 1 (in  MO frame)
    r_f_probe1     = r_measure+dr_p1_MO;                         % [m] farthest clock distance to MO
    r_n_probe1     = r_measure;                                  % [m] nearest clock distance to MO
    dtn_dtf_probe1 = frames_dtn_dtf(M_MO,r_f_probe1,r_n_probe1); % [-] clock differential
    g_m_probe1     = gravimeter(dtn_dtf_probe1,dr_orb_MO);       % [m/s^2] measured g
    
    % determine effects on gravimeter from probe 2 (in  MO frame)
    r_f_probe2     = r_measure+dr_p2_MO;                         % [m] farthest clock distance to MO
    r_n_probe2     = r_measure;                                  % [m] nearest clock distance to MO
    dtn_dtf_probe2 = frames_dtn_dtf(M_MO,r_f_probe2,r_n_probe2); % [-] clock differential
    g_m_probe2     = gravimeter(dtn_dtf_probe2,dr_orb_MO);       % [m/s^2] measured g
    
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
plot([v_obj_all(1) v_obj_all(end)], [r_2_gravObj(M_MO,r_measure) r_2_gravObj(M_MO,r_measure)],'k--','LineWidth',2)

% clean up plot
legend('Orbital Gravitmeter','USF Faster Gravimeter','USF Slower Gravimeter','Truth Gravity','FontSize',16,'location','NW');
xlabel({'Dimensional Velocity of Massed Object [fraction of c]','With Respect to USF'},'FontSize',16);
ylabel('Measured $g(r)~\left[\frac{m}{s^2}\right]$','FontSize',16,'Interpreter','latex');
grid on
xticks([0:.1:1]);
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'fontsize',16)
annotation(fig, 'textbox', [.13 .10 .8 .2], 'String'...
    ,sprintf('Mass of massed object: %d [Solar Masses]',M_MO/Ms)...
    ,'EdgeColor','none','FontSize',14);
annotation(fig, 'textbox', [.13 .07 .8 .2], 'String'...
    ,sprintf('Gravimeter clocks original distance apart: %d [km]',dr_orb_MO_0/1e3)...
    ,'EdgeColor','none','FontSize',14);
annotation(fig, 'textbox', [.13 .04 .8 .2], 'String'...
    ,sprintf('Measurement distance from center of mass: %0.1f [AU]',r_measure/AU)...
    ,'EdgeColor','none','FontSize',14);
annotation(fig, 'textbox', [.13 .01 .8 .2], 'String'...
    ,sprintf('Speed of probes relative to massed object: %0.1f [fraction of c]',probe_dv)...
    ,'EdgeColor','none','FontSize',14);

%% supporting function
    function dtn_dtf = frames_dtn_dtf(M,r_f,r_n)
        % (in MO frame)
        dt_f      = gamma_inv_P(M,r_f);    % time dilation of clock farthest from MO
        dt_n      = gamma_inv_P(M,r_n);    % time dilation of clock nearest to MO
        dtn_dtf   = dt_n/(dt_f);    % relative time differential between closest and farthest clock
    end
    
    function r2 = solve_for_r2(M,r1,dr)
        % initial guess
        r2_upper = r1 + 2*dr;
        r2_lower = r1;
        r2 = (r2_lower + r2_upper)/2;
        dr_guess = prop_dist(M,r1,r2);
        error = dr_guess - dr;
        while (1e-9 < abs(error) || dr/(2^25) > abs(r2_upper-r2_lower))
            if 0 < error
                r2_upper = r2;
            else
                r2_lower = r2;
            end
            r2 = (r2_lower + r2_upper)/2;
            dr_guess = prop_dist(M,r1,r2);
            error = dr_guess - dr;
        end
    end
end