%% clear close
clear
close all

%% constants
c =  299792458; % [m/s] speed of light
h = 6.62607015e-34;% [m^2 kg / s] Planck's Constant
kB = 5.670374419184429453970e-8; % [J?m?2?s?1?K?4] stefan-boltzmann constant
b = 2.897771955e-3; % [m?K] Wien's displacement constant
G = 6.6744e-11; % [m^3/(kg s)] gravitational constant
re = 6371000; % [m] earth's mean radius
Me = 5.97219e24; % [kg] earth's mass
Ms = 333000*Me; % [kg] sun's mass
eMax = c^2/2;
vGalaxy = 0.581152e6; %[m/s] how fast our galaxy is moving
planckl = @(lambda,T) 2*h*c^2./(lambda.^5.*(exp(h*c./(kB*lambda.*T))-1));
planckl_odd = @(lambda,T) 2*h*c^2./(lambda.^2.*(exp(h*c./(kB*lambda.*T))-1));
planckf = @(v,T) 2*h*v.^3./(c^2*(exp(h*v/(T))-1));
ly_per_parsec = 26/8;
ls_per_ly = 365*24*60*60;
m_per_ls = c;
m_per_parsec =  (m_per_ls * ls_per_ly * ly_per_parsec);
c_parsec = c / m_per_parsec;
G_parsec = G / (m_per_ls * ls_per_ly * ly_per_parsec)^3;

v_green = .55e-6;
m_green = v_green * h;
m_electron = 1e-31;
photon_per_electron = m_electron / m_green;

%% taylor series expantion of 1/gamma error
v_c = [0:1e-3:1];
gamma = sqrt(1-v_c.^2);
gamma_t1 = ones(size(v_c));
gamma_t2 = 1./(ones(size(v_c)) + 0.5*v_c.^2);
gamma_t3 = ones(size(v_c));
gamma_t6 = ones(size(v_c));
n = -.5;
p = n;
k = 100;
num = n;
for i = 1 : k - 1
    fac_i = factorial(i);
    c_i = num/fac_i;
    gamma_t6 = gamma_t6 + abs(c_i)*v_c.^(2*i);
    if i <=2
        gamma_t3 = gamma_t3 + abs(c_i)*v_c.^(2*i);
    end
    p = n - i;
    num = num*p;
end
gamma_t3 = 1./gamma_t3;
gamma_t6 = 1./gamma_t6;
figure()
subplot(2,1,1)
plot(v_c,gamma,'k','linewidth',2);
hold on
plot(v_c,gamma_t1,'-r','linewidth',2);
plot(v_c,gamma_t2,'--r','linewidth',2);
plot(v_c,gamma_t3,'-.r','linewidth',2);
plot(v_c,gamma_t6,':r','linewidth',2);
plot([0 1],[gamma_t6(end) gamma_t6(end)],'--k','LineWidth',2)
legend('$\frac{1}{\gamma}$','$1^{st}~Order$','$2^{nd}~Order$','$3^{rd}~Order$'...
    ,'$1000^{th}~Order$','$\frac{1}{\gamma}$ Estimate @ $\frac{1}{\gamma}=1$','Location','SW'...
    ,'fontsize',14,'interpreter','latex')
xlabel('v/c','FontSize',20,'interpreter','latex');
ylabel('$\frac{1}{\gamma}$ Estimate','FontSize',20,'interpreter','latex');
title({'$\frac{1}{\gamma}$ Estimates','Taylor Series Evaluated at v/c = 0'}...
    ,'FontSize',16,'interpreter','latex')
grid on
xticks([0:.1:1])
xticklabels([0:.1:1])
yticks([0:.1:1])
yticklabels([0:.1:1])
subplot(2,1,2)
plot(v_c,(gamma_t1-gamma),'-r','linewidth',2);
hold on
plot(v_c,(gamma_t2-gamma),'--r','linewidth',2);
plot(v_c,(gamma_t3-gamma),'-.r','linewidth',2);
plot(v_c,(gamma_t6-gamma),':r','linewidth',2);
plot([0 1],[gamma_t6(end)-gamma(end) gamma_t6(end)-gamma(end)],'--k','LineWidth',2)
legend('$1^{st}~Order$','$2^{nd}~Order$','$3^{rd}~Order$','$1000^{th}~Order$'...
    ,'$\frac{1}{\gamma}$ Error @ $\frac{1}{\gamma}=1$','Location','NW','fontsize',14,'interpreter','latex')
xlabel('v/c','FontSize',20,'interpreter','latex');
ylabel('$\frac{1}{\gamma}$ Error','FontSize',20,'interpreter','latex');
title({'$\frac{1}{\gamma}$ Error','Taylor Series Evaluated at v/c = 0'}...
    ,'FontSize',16,'interpreter','latex')
grid on
xticks([0:.1:1])
xticklabels([0:.1:1])
yticks([0:.1:1])
yticklabels([0:.1:1])

%% graph twin ages more with same acceleration when there is more travel
vfc = 0.5;
ax = 0.5*(vfc)^2;
gamma_inv = sqrt(1-ax/.5);
tp = [0:1e-2:3]*365*24*60*60;
x = vfc.*tp;
t = gamma_inv.*tp;
delta_t = tp - t;
fig = figure(1)
subplot(2,1,1);
plot(x/(365*24*60*60),delta_t/(365*24*60*60),'k','LineWidth',2);
xlabel('Distance Traveled (Light-Years)','FontSize',20);
ylabel({'Age Difference','(Years)'},'FontSize',20);
xticks([0:.1:3*vfc]);
yticks([0:.1:3]);
title('Traveling Twin at $\frac{1}{2}c$','Interpreter','latex','FontSize',20)
grid on
subplot(2,1,2);
plot(zeros(size(delta_t)),delta_t/(365*24*60*60),'k','LineWidth',2);
xlabel('Change in Acceleration From Original Profile','FontSize',20);
ylabel({'Age Difference','(Years)'},'FontSize',20);
yticks([0:.1:3]);
grid on
close(fig);

%% velocity is not the cause of time dilation
v1 = [0:.01:1];
v2 = [v1];
vr = (v1+v2)./(1+(v1.*v2));
fig = figure(1);
plot(vr,zeros(size(vr)),'k','LineWidth',2);
xlabel('Relative Velocity (Fraction of $c$)','FontSize',20,'Interpreter','latex');
ylabel('$\frac{dt_2}{dt_3}$','Interpreter','latex','FontSize',20);
xticks([0:.1:3*vfc]);
grid on


%% refraction correction
v = c*[0:1e-3:1-1e-3];
n_orig = c./v;
gamma_inv = sqrt(1-(v/c).^2);
n_new = flip(1./gamma_inv);
plot(n_orig, n_new);

%% circular orbit speed with relativity
Mg = 1e11*Ms; % pass of center of galaxy
rs  = G_parsec*Mg / c_parsec^2;
r = rs:1e-6:10;
Ve_c = m_per_parsec*sqrt(G_parsec*Mg ./ r)/c;
dt_dtp_P = sqrt(1-rs./r);
Ve_c_observed = Ve_c.*dt_dtp_P;
plot(r,Ve_c_observed);
set(gca,'XScale','log');

%% light-g ratio for different things
photons_Sun_per_Sec = 1e45;
obj.E.T = 287.039; % [K] average earth temperature
% avg temps
obj.sun.T = 5778;  % [K] average sun temperature
obj.V.T = 737;     % [K] average venus temperature
obj.E.T = 287.039; % [K] average earth temperature
obj.M.T = 210.372; % [K] average mars temperature
% https://www.google.com/url?sa=t&rct=j&q=&esrc=s&source=web&cd=&cad=rja&uact=8&ved=2ahUKEwi4hdjByaj5AhXBGVkFHUFJAgUQFnoECBMQAw&url=https%3A%2F%2Fwww.smartconversion.com%2F(X(1))%2FotherInfo%2FTemperature_of_planets_and_the_Sun.aspx&usg=AOvVaw0FW8wH3CVfbSVkNm2ae-Z2
obj.J.T = 165;     % [K] average Jupiter temperature
obj.sat.T = 134;   % [K] average Saturn temperature
obj.U.T = 76;      % [K] average Uranus temperature

% avg radius
obj.sun.r = 696.34e6; % [m] average sun radius 
obj.V.r = 6.0518e6;   % [m] average venus radius 
obj.E.r = re;         % [m] average earth radius 
obj.M.r = 3.3895e6;   % [K] average mars radius
obj.J.r = 69.911e6;   % [K] average Jupiter radius
obj.sat.r = 58.232e6; % [K] average Saturn radius
obj.U.r = 25.362e6;   % [K] average Uranus radius

% avg mass
obj.sun.M = 1.989e30; % [m] average sun radius 
obj.V.M = 4.867e24;   % [m] average venus radius 
obj.E.M = Me;         % [m] average earth radius 
obj.M.M = 6.39e23;   % [K] average mars radius
obj.J.M = 1.898e27;   % [K] average Jupiter radius
obj.sat.M = 5.683e26; % [K] average Saturn radius
obj.U.M = 8.681e25;   % [K] average Uranus radius

% get peak freq, power, power density, g, pwr den/g ratio
obj_names = fieldnames(obj);
figure(1);
hold off
colors = {[0 0 0]...
    ,     [1 1 0]...
    ,     [1 0 1]...
    ,     [0 1 1]...
    ,     [1 0 0]...
    ,     [0 1 0]...
    ,     [0 0 1]...
    };
for io = 1 : length(obj_names)
    obj_name = obj_names{io};
    lamb_peak = (b / obj.(obj_name).T);
    obj.(obj_name).f = c / (b / obj.(obj_name).T); % weins dispacement law to get peak freq
    obj.(obj_name).P  = 4*pi*obj.(obj_name).r^2 * obj.(obj_name).T^4 * kB;%integral(@(x) planckf(x,obj.(obj_name).T),obj.(obj_name).f/100,obj.(obj_name).f*100);%
    r_range = [obj.(obj_name).r:1000:obj.(obj_name).r+1e7];
    obj.(obj_name).Pd = obj.(obj_name).P ./ (4*pi*r_range.^2);
    obj.(obj_name).g  = G * obj.(obj_name).M ./ r_range.^2;
    obj.(obj_name).g2pd = obj.(obj_name).g ./ obj.(obj_name).Pd;
    obj.(obj_name).pd2g2f = (obj.(obj_name).Pd/obj.(obj_name).T^3) ./ obj.(obj_name).g;
    obj.(obj_name).EK_avg = obj.(obj_name).T*kB*3/2;
    obj.(obj_name).eK_avg = obj.(obj_name).T*kB*3/(2*obj.(obj_name).M);
    yyaxis right;
    plot(r_range-obj.(obj_name).r,obj.(obj_name).g2pd*1e3,'-','Color',colors{io},'linewidth',1);
    yyaxis left;
    plot(0,obj.(obj_name).T^3/1e8,'x','Color',colors{io},'linewidth',1);
    hold on
end
legend(obj_names,'location','best');
yyaxis right;
ylabel('Power Density to Gravitational Acceleration Ratio');
yyaxis left;
ylabel('Temperature [K]');
xlabel('Distant (r) away from Surface');
set(gca,'YScale','log');
yyaxis right;
set(gca,'YScale','log');
%% wave path length for different light wave lengths
f_green_0 = 5.45e15;
lambda_green_0 = c/f_green_0;
t0 = 0;
tend = 1 / f_green_0;
t = t0:tend/1e3:tend;
wave0 = @(t) sin(2*pi*f_green_0*t);
signal0 = wave0(t);
[~,idx_max] = max(signal0);
lambda_est = c*t(idx_max)*4;
delta_v = c/5;
dx_dxp=sqrt(1-(delta_v/c)^2);
tau = delta_v/c;
n = 1/(1-tau);
v0 = c;
v1 = v0/n;
lambda_green_1 = v1 / f_green_0;
f_green_1 = c/lambda_green_1;
wave1 = @(t) sin(2*pi*f_green_0*t*n);
signal1 = wave1(t);

figure(1);
plot(t*c,signal0);
hold on
plot(t*c,signal1,'--');

S0 = sum(sqrt((signal0(2:end) - signal0(1:end-1)).^2 + (c*(t(2)-t(1)))^2));
S1 = sum(sqrt((signal1(2:end) - signal1(1:end-1)).^2 + (c*(t(2)-t(1)))^2));

%% test different dx_dxp functionsx
dt_dtp = [0.1:1e-6:1];
v1_p = sqrt(1-dt_dtp.^2);
dx_dxp = sqrt(1-v1_p.^2);
dx_dxp_weird = (1-v1_p.^2);
dx_dt = dx_dxp_weird ./ dt_dtp;

figure(1);
subplot(3,1,1)
plot(v1_p,v1_p.*dx_dxp_weird./dt_dtp,'k','LineWidth',3);
grid on
xlabel('$\frac{v''}{c}$','interpreter','latex','fontsize',20);
ylabel('$\frac{v}{c}$','interpreter','latex','fontsize',20);
title('$\frac{dx}{dx''}<\frac{1}{\gamma}$','interpreter','latex','fontsize',20);
xticks([0:.1:1]);
yticks([0:.1:1]);

subplot(3,1,2)
plot(v1_p,v1_p.*dx_dxp./dt_dtp,'k','LineWidth',3);
grid on
xlabel('$\frac{v''}{c}$','interpreter','latex','fontsize',20);
ylabel('$\frac{v}{c}$','interpreter','latex','fontsize',20);
title('$\frac{dx}{dx''}=\frac{1}{\gamma}$','interpreter','latex','fontsize',20);
xticks([0:.1:1]);
yticks([0:.1:1]);

subplot(3,1,3)
plot(v1_p,v1_p./dt_dtp,'k','LineWidth',3);
grid on
xlabel('$\frac{v''}{c}$','interpreter','latex','fontsize',20);
ylabel('$\frac{v}{c}$','interpreter','latex','fontsize',20);
title('$\frac{dx}{dx''}>\frac{1}{\gamma}$','interpreter','latex','fontsize',20);
hold on
plot([0 1],[1 1],'r--','LineWidth',2)
xticks([0:.1:1]);
yticks([0:1:10]);
theselabels = yticklabels;
theselabels{2} = 'c';
yticklabels(theselabels);

%% alpha centauri trip distance changed
vp = 0;
dt_dtp = 0.5;
v1_p = sqrt(1-dt_dtp^2);
distAS_p = 4; % [light years] approx distance to alpha centauri
time_p = distAS/v1_p; % [years] time it takes to get there
time_m = time_p * dt_dtp;
distAS_m = time_m * v1_p;
dx_dxp_mesured = distAS_m / distAS_p;
delta_error = dx_dxp_mesured - dt_dtp;
delta_vel = v1_p - V1_m_if;

%% GPS test
dtGPS = 1; % GPS dt'
GPSfaster = 45e-6/(24*60^2); % 45us/day faster
dte = dtGPS - dtGPS*GPSfaster; % earth dt
deltaGPS_t = dte-dtGPS; % 45us/day slower
delta_GPS_r = 1e3; % [m] GPS altitude above ground level
rGPS = re + delta_GPS_r; % [m] GPS distance from center of earth
drGPS = delta_GPS_r/1e5;
r = re:drGPS:rGPS;
g = Me*G./r.^2;
dSE = g.*drGPS;
deltaSE = cumsum(dSE);
dSEend = deltaSE(end);
dt1_dt2 = dtGPS*sqrt(1-dSEend/eMax); % agrees with deltaGPS_t in line 13
grad_t = (dt1_dt2-dtGPS)/delta_GPS_r;
dv = delta_GPS_r/dtGPS;
a_est = (eMax/delta_GPS_r)*(1-(dv*grad_t+1)^2);
gGPS = G*Me/rGPS^2;
ge = G*Me/re^2;
g_est = G*Me*(1/re - 1/rGPS)/delta_GPS_r;%(gGPS^dtGPS*ge^dte)^(1/(dtGPS + dte)); % matches a_est line 25
d_error = g_est - a_est;
p_error = 100*d_error/g_est;

%% mean gradient test
fr = @(x) 1 ./ x.^2;
mean_fr = integral(fr,2,3);
geoMean_fr = sqrt(fr(2)*fr(3));
%% twins travel together first
dist = 1; %lt-sec
v0 = 0.2; % percent c for both twins at start wrt original reference frame
v1 = 0; % "moving" twin drops back to 0 wrt original reference frame (-0.5 wrt "stationary" twin)
v2_0 = 0.2; % how fast (percent c) the "moving" twin catches back up from the twins' perspective
v2 = (v0 + v2_0) / (1 + v0*v2_0); % how fast (percent c) the "moving" twin catches back up from the original frame's perspective

dts_dtp = sqrt(1-.5); % "stationary" twin's time dilation wrt original reference frame
dt1_dtp = sqrt(1-v1); % "moving" twin's time dilation wrt original reference frame at t1
dt2_dtp = sqrt(1-v2); % "moving" twin's time dilation wrt original reference frame at t2
t_away = dist / v0; % how much time passes for away trip (t1-t0), in original reference frame time, to cover distance
t_return = dist / (v2 - v0);  % how much time passes for return trip (t2-t1), in original reference frame time, covers original distance and then some because "stationary" twin is moving at 0.5c

total_dts = sum(dts_dtp*[t_away  t_return]); % total passing of time, from start of twin's separation, for "stationary" twin, as measured by "stationary" twin
total_dtm = sum([dt1_dtp*(t_away)  dt2_dtp*(t_return)]); % total passing of time, from start of twin's separation, for "moving" twin, as measured by "moving" twin
ts = (t_away + t_return)*(1+v0*v0)*dts_dtp^-1;
tm = t_return*(1+v2*v2)*dt2_dtp^-1;

%% space differential
v = 0.5; % fraction of c
fc = 1; % fraction of c

dt_dtp = sqrt(1-v);
time = 1; % light-second
dp = 1;
d = (fc-v) * time;
dx_dxp = (d / dt_dtp);

%% space differential 2
v = 2; % fraction speed of light
gamma = 1/sqrt(1-v^2);
dx_dxp = (fc/(v*gamma));
dx_dxp2 = sqrt(fc^2/v^2-1);

%% testing solar gravity as light
a_e2s =  0.006; % [m/s^2] acceleration of earth 2 sun
dr  = 1e3; % [m] precision
dt  = sqrt(2*dr/a_e2s);
SW = a_e2s*dr;
dt_dtp = sqrt(1 - SW / c^2);

rSunE = 152.03e9; % [m] distance from sun to earth
de2 = rSunE + dr; % [m] distance from sun to earth plus precision
P_sun = 3.9e26; % [W] [kg m^2/s^2 /s] energy during duration
E_sun = P_sun * dt; % [kg m^2/s^2] energy during duration
Ade1 = 4*pi*rSunE^2;
Ade2 = 4*pi*de2^2;
E_e1 = E_sun / (Ade1);
E_e2 = E_sun / (Ade2);
dE_e = E_e1 - E_e2;
a_e2s_est = dE_e / dr;

%% solar time derivative gradient
dtSunE_far = 1;
delta_SunE_r = eMax; % [m] Additional distance
rSunE_far = rSunE + delta_SunE_r; % [m] earth distance from center of Sun
drSunE = delta_SunE_r/c;
r = rSunE:drSunE:rSunE_far;
g = Ms*G./r.^2;
dSE = g.*drSunE;
deltaSE = cumsum(dSE);
dSEend = deltaSE(end);
dt1_dt2 = dtSunE_far*sqrt(1-dSEend/eMax); % agrees with deltaGPS_t in line 13
grad_t = (dt1_dt2-dtSunE_far)/delta_SunE_r;
dv = delta_SunE_r/dtSunE_far;
a_est = (eMax/delta_SunE_r)*(1-(dv*grad_t+1)^2);
gSunE_far = G*Ms/rSunE_far^2;
gSunE = G*Ms/rSunE^2;
g_est = sqrt(gSunE_far*gSunE);%(gGPS^dtGPS*ge^dte)^(1/(dtGPS + dte)); % matches a_est line 25
d_error = g_est - a_est;
p_error = 100*d_error/g_est;

%% total specific energy time dilation
ep_emax = 1/3;
ek_emax = sqrt(2);
gamma_inv_SQ_ep = 1 - ep_emax;
gamma_inv_SQ_ek = 1 - ek_emax;
tau_SQ_ep = 1 - gamma_inv_SQ_ep;
tau_SQ_ek = 1 - gamma_inv_SQ_ek;
tau_SQ_eT = tau_SQ_ep + tau_SQ_ek;
gamma_inv_SQ_eT = 1 - tau_SQ_eT;
gamma_inv_SQ_eT2 = 1 - ep_emax - ek_emax;
delta = gamma_inv_SQ_eT - gamma_inv_SQ_eT2;

%% speed changes imperceptably
v = [0.9:1e-5:1-1e-5];
vp = 2*v./(1+v.^2);
plot(v,vp,'k','LineWidth',3)
xlabel('Velocity of Objects Going in Opposite Directions','FontSize',16);
ylabel('Velocity of Moving Objects as Seen By Moving Objects','FontSize',16);
grid on
title('Changes In Speed Becomes Imperceptable When Approaching c','FontSize',16)

%% mass or speed of a photon
lambda1 = 1e-6;
lambda2 = 11e-6;
E1 = h*c/lambda1;
E2 = h*c/lambda2;
m1 = 2 * h / (lambda1*c);
c2  = 2 * h / (lambda2*m1);
m2 = 2 * h / (lambda2*c);

%% GPS


disp("done");
