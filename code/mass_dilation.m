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
gd1     = 2*G*m/(d);  % [J/kg]   initial relative xspecific potential energy
v       = 0.6;        % [-]      fraction of the speed of light of orbs
gamma_v = gamma(v);   % [-]      1/sqrt(1-v^2/c^2)

% initialize other variables
dy  = (d-d_min)/1e4;     % increment steps to numerical solution
ds  = d:-dy:d_min;       % all numerical steps
gds = ones(size(ds))*gd1;% specific potential energy
vs  = zeros(size(ds));   % relative velocity of orbs
ts  = zeros(size(ds));   % proper time passed

% incremental solution of orb pairs relative velocity and time passed
for id = 2 : length(ds)
    % this relative specific potential energy for orbs
    gds(id) = 2*G*m/(ds(id));
    
    % delta relative specific potential energy for orbs
    delta_gd = gds(id)-gd1;
    
    % relative velocity between them
    vs(id) = sqrt(2*delta_gd);
    
    % time for distance to close by mean relative velocity
    ts(id) = ts(id-1) + dy/mean([vs(id),vs(id-1)]);
end

% total passage of proper time until orbs contact in years and months
total_time_months = max(ts)*seconds2months;

%% Stationary orbs
my      = m/gamma_v^2; % [kg]   mass of stationary orb is traveling orb's mass
gd1_my  = 2*G*my/(d);  % [J/kg] initial specific potential energy

% time passed, as measured by stationary orbs
ts_gamma = ts*gamma_v;

% total passage of proper time until orbs contact in years and months
total_time_months_my  = max(ts_gamma)*seconds2months;

% initialize stationary orbs with mass my distance steps
dy_my = dy;             % increment steps to numerical solution
ds_my = d:-dy_my:d_min; % all numerical steps

% initialize other variables
vs_my  = zeros(size(ds_my));
gds_my = ones(size(ds_my))*gd1_my;
ts_my  = zeros(size(ds_my));

% incremental solution of orb pairs relative velocity and time passed
for id = 2 : length(ds_my)
    % this relative specific potential energy
    gds_my(id) = 2*G*my/(ds_my(id));
    
    % delta relative specific potential energy
    delta_gd_my = gds_my(id)-gd1_my;
    
    % relative velocity between them
    vs_my(id) = sqrt(2*delta_gd_my);
    
    % time for distance to close by mean relative velocity
    ts_my(id) = ts_my(id-1) + dy_my/mean([vs_my(id),vs_my(id-1)]);
end

%% Plot Results
figure(1);
% plot the movement of each orb makes towards its pair
subplot(2,1,1)
plot((d-ds)/2,(d-interp1(ts_my,ds_my,ts_gamma))/2,'-b','LineWidth',1.5)
xlim([0 d/2]);
ylim([0 d/2]);
grid on
xlabel('Stationary Orbs','FontSize',20);
ylabel('Traveling Orbs','FontSize',20);
title({'Distance Each Orb Traveled Towards The Other [m]'},'fontsize',16);

% plot the percent difference in movement between pairs of orbs
percent_difference = 100*abs((d-interp1(ts_my,ds_my,ts_gamma))/2 - (d-ds)/2)./(dy);
subplot(2,1,2)
plot((d-ds)/2,percent_difference,'-b','LineWidth',2)
xlim([0 d/2]);
ylim([0 100]);
grid on
xlabel({'Distance Each Orb in Stationary Frame'...
    ,'Traveled Towards The Other [m]'},'FontSize',20);
ylabel({'\%$\Delta$'},'FontSize',20,'Interpreter','latex');
title({'\%$\Delta=100\times\frac{|Traveling-Stationary|}{precision}$'}...
    ,'Interpreter','latex','fontsize',16);

% print ellapsed proper (AAK wall) time for each pair or orbs
fprintf('Elapsed Time for Traveling Orbs: %0.1f [months]\n',total_time_months);
fprintf('Elapsed Time for Stationary Orbs: %0.1f [months]\n',total_time_months_my);