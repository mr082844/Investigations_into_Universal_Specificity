%% constants and functions
G             = 6.6744e-11; % [m^3/(kg s)] gravitational constant
gamma         = @(v) 1./sqrt(1-v.^2);
seconds2years = 1/60^2/24/365;
years2months  = 12;

%% Traveling orbs
% initial conditions
rho     = 22000;          % [kg/m^3] density of osmium
r       = 1e-1;           % [m]      radius of each orb
vol     = 4*pi*r^3/3;     % [m^3]    volume of each orb
m       = rho*vol;        % [kg]     mass of each orb
d       = 1e2;            % [m]      initial distance between orbs' surfaces
g1      = 2*G*m/(d+2*r);  % initial relative acceleration for orbs
v       = 0.6;            % fraction of the speed of light of orbs
gamma_v = gamma(v);       % 1/sqrt(1-v^2/c^2)

% initialize other variables
dy = d/1e4;             % increment steps to numerical solution
ds = d:-dy:0;           % all numerical steps
gs = ones(size(ds))*g1; % relative acceleration
vs = zeros(size(ds));   % relative velocity of orbs
ts = zeros(size(ds));   % proper time passed

% incremental solution of orb pairs relative velocity and time passed
for id = 2 : length(ds)
    % this relative acceleration for orbs
    gs(id) = 2*G*m/(ds(id)+2*r);
    
    % delta relative acceleration for orbs
    dg = gs(id)-g1;
    
    % relative velocity between them
    vs(id) = sqrt(2*dg);
    
    % time for distance to close by mean relative velocity
    ts(id) = ts(id-1) + dy/mean([vs(id),vs(id-1)]);
end

% total passage of proper time until orbs contact in years and months
total_time_years = max(ts)*seconds2years;
total_time_months = total_time_years*years2months;

%% Stationary orbs
m0      = m/gamma_v^2;    % mass of stationary orb is predicted rest mass of traveling orbs
g1_m0   = 2*G*m0/(d+2*r); % initial relative acceleration for stationary orbs

% time passed, as measured by stationary orbs
ts_gamma = ts*gamma_v;

% total passage of proper time until orbs contact in years and months
total_time_years_m0  = max(ts_gamma)*seconds2years;
total_time_months_m0 = total_time_years_m0*years2months;

% initialize stationary orbs with mass m0 distance steps
dy_m0 = dy;         % increment steps to numerical solution
ds_m0 = d:-dy_m0:0; % all numerical steps

% initialize other variables
vs_m0 = zeros(size(ds_m0));
gs_m0 = ones(size(ds_m0))*g1_m0;
ts_m0 = zeros(size(ds_m0));

% incremental solution of orb pairs relative velocity and time passed
for id = 2 : length(ds_m0)
    % this relative acceleration
    gs_m0(id) = 2*G*m0/(ds_m0(id)+2*r);
    
    % delta relative acceleration
    dg_m0 = gs_m0(id)-g1_m0;
    
    % relative velocity between them
    vs_m0(id) = sqrt(2*dg_m0);
    
    % time for distance to close by mean relative velocity
    ts_m0(id) = ts_m0(id-1) + dy_m0/mean([vs_m0(id),vs_m0(id-1)]);
end

%% Plot Results
% plot the movement of each pair of orbs with respect to each other
figure(1);
subplot(2,1,1)
hold off
plot(ds/2,interp1(ts_m0,ds_m0,ts_gamma)/2,'-b','LineWidth',1.5)
title({'Distance Orbs Traveled Towards Each Other [m]'},'fontsize',16);
grid on
xlabel('Stationary Orbs','FontSize',20);
ylabel('Traveling Orbs','FontSize',20);

% plot the percent difference in movement between pairs of orbs
percent_difference = abs(interp1(ts_m0,ds_m0,ts_gamma) - ds)./(dy);
subplot(2,1,2)
hold off
plot(ds/2,100*percent_difference,'-b','LineWidth',2)
title({'% Difference'},'fontsize',16);
grid on
xlabel('Stationary Orbs'' Distance Moved [m]','FontSize',20);
ylabel({'$100\times\frac{|Traveling-Stationary|}{precision}$'},'FontSize',20,'Interpreter','latex');
ylim([0 100]);

disp("Done.");