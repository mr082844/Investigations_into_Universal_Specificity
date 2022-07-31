%% constants
c =  299792458; % [m/s] speed of light
h = 6.62607015e-34;% [m^2 kg / s] Planck's Constant
G = 6.6744e-11; % [m^3/(kg s)] gravitational constant
re = 6371000; % [m] earth's mean radius
Me = 5.97219e24; % [kg] earth's mass
Ms = 333000*Me; % [kg] sun's mass
eMax = c^2/2;

%% GPS test
dtGPS = 1; % GPS dt'
GPSfaster = 45e-6/(24*60^2); % 45us/day faster
dte = dtGPS - dtGPS*GPSfaster; % earth dt
deltaGPS_t = dte-dtGPS; % 45us/day slower
delta_GPS_r = 1e1; % [m] GPS altitude above ground level
rGPS = re + delta_GPS_r; % [m] GPS distance from center of earth
drGPS = delta_GPS_r/1000;
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
v0 = 0.9; % percent c for both twins at start wrt original reference frame
v1 = .7; % "moving" twin drops back to 0 wrt original reference frame (-0.5 wrt "stationary" twin)
v2_0 = 0.4; % how fast (percent c) the "moving" twin catches back up from the twins' perspective
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
