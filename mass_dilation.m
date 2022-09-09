%% constants and functions
G               = 6.6744e-11; % [m^3/(kg s)] gravitational constant
gamma           = @(v) 1./sqrt(1-v.^2);
seconds2months  = 12/60^2/24/365;

%% Traveling orbs
% initial conditions
rho     = 22000;      % [kg/m^3] density of osmium
r       = 1e-1;       % [m]      radius of each orb
vol     = 4*pi*r^3/3; % [m^3]    volume of each orb
m       = rho*vol;    % [kg]     mass of each orb
d       = 1e2;        % [m]      initial distance between orbs' surfaces
d_min   = 2*r;        % [m]      minimum distance between center mass of orbs
gd1     = 2*G*m/(d);  % [J/kg]   initial specific potential energy
v       = 0.6;        % [-]      fraction of the speed of light of orbs
gamma_v = gamma(v);   % [-]      1/sqrt(1-v^2/c^2)

% initialize other variables
dy  = (d-d_min)/1e4;     % increment steps to numerical solution
ds  = d:-dy:d_min;       % all numerical steps
gds = ones(size(ds))*gd1; % relative acceleration
vs  = zeros(size(ds));   % relative velocity of orbs
ts  = zeros(size(ds));   % proper time passed

% incremental solution of orb pairs relative velocity and time passed
for id = 2 : length(ds)
    % this relative acceleration for orbs
    gds(id) = 2*G*m/(ds(id));
    
    % delta relative acceleration for orbs
    delta_gd = gds(id)-gd1;
    
    % relative velocity between them
    vs(id) = sqrt(2*delta_gd);
    
    % time for distance to close by mean relative velocity
    ts(id) = ts(id-1) + dy/mean([vs(id),vs(id-1)]);
end

% total passage of proper time until orbs contact in years and months
total_time_months = max(ts)*seconds2months;

%% Stationary orbs
m0      = m/gamma_v^2; % [kg]   mass of stationary orb is traveling orb's rest mass
gd1_m0  = 2*G*m0/(d);  % [J/kg] initial specific potential energy

% time passed, as measured by stationary orbs
ts_gamma = ts*gamma_v;

% total passage of proper time until orbs contact in years and months
total_time_months_m0  = max(ts_gamma)*seconds2months;

% initialize stationary orbs with mass m0 distance steps
dy_m0 = dy;             % increment steps to numerical solution
ds_m0 = d:-dy_m0:d_min; % all numerical steps

% initialize other variables
vs_m0  = zeros(size(ds_m0));
gds_m0 = ones(size(ds_m0))*gd1_m0;
ts_m0  = zeros(size(ds_m0));

% incremental solution of orb pairs relative velocity and time passed
for id = 2 : length(ds_m0)
    % this relative acceleration
    gs_m0(id) = 2*G*m0/(ds_m0(id));
    
    % delta relative acceleration
    delta_gd_m0 = gs_m0(id)-gd1_m0;
    
    % relative velocity between them
    vs_m0(id) = sqrt(2*delta_gd_m0);
    
    % time for distance to close by mean relative velocity
    ts_m0(id) = ts_m0(id-1) + dy_m0/mean([vs_m0(id),vs_m0(id-1)]);
end

%% Plot Results
figure(1);
% plot the movement of each orb makes towards its pair
subplot(2,1,1)
hold off
plot((d-ds)/2,(d-interp1(ts_m0,ds_m0,ts_gamma))/2,'-b','LineWidth',1.5)
title({'Distance Each Orb Traveled Towards The Other [m]'},'fontsize',16);
grid on
xlabel('Stationary Orbs','FontSize',20);
ylabel('Traveling Orbs','FontSize',20);

% plot the percent difference in movement between pairs of orbs
percent_difference = 100*abs((d-interp1(ts_m0,ds_m0,ts_gamma))/2 - (d-ds)/2)./(dy);
subplot(2,1,2)
hold off
plot((d-ds)/2,percent_difference,'-b','LineWidth',2)
title({'% Difference'},'fontsize',16);
grid on
xlabel({'Distance Each Orb in Stationary Frame'...
    ,'Traveled Towards The Other [m]'},'FontSize',20);
ylabel({'$100\times\frac{|Traveling-Stationary|}{precision}$'}...
    ,'FontSize',20,'Interpreter','latex');
ylim([0 100]);

% print ellapsed proper (AAK wall) time for each pair or orbs
fprintf('Elapsed Time for Traveling Orbs: %0.1f [months]\n',total_time_months);
fprintf('Elapsed Time for Stationary Orbs: %0.1f [months]\n',total_time_months_m0);