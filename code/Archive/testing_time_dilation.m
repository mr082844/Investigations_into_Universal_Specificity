r1 = 50;
r2 = 25;
M = 5.97219e24; % kg earth
G = 6.6743e-11; % m3 kg-1 s-2
c = 299792458; % m / s
rearth = 6.371e6; % m mean earth radius
ralpha = 4.132e16; % m to alpha centauri from earth
ralpha = 4.132e16; % order of magnitude less

% bunch of rs for position 1 and 2
rs1 = [80:1:500];
rs2 = rs1-40;
dr = 1e9;
rs_inf = [rearth:dr:ralpha];

SPE = -dr*G*M./(rs_inf.^2);
cus_SPE = cumsum(SPE);
dt_dtp = sqrt(1-cus_SPE/(.5*c^2));
plot(rs_inf,1./dt_dtp);

% time dilation 1
dt1_dtp = sqrt(1-2*G*M./(rs1*c^2));
rs1p = rs1 ./ dt1_dtp;

% time dilation 2 and apparent length contraction
dt2_dtp = sqrt(1-2*G*M./(rs2*c^2));
rs2p = rs2 ./ dt1_dtp;

% time dilation 2 / 1
Gp = G ./ dt1_dtp.^2;
dt2_dt1 = dt2_dtp./dt1_dtp;
dt2_dt1p = sqrt(1-2*Gp*M./((rs1p - rs2p)*c^2));

figure(1)
hold off
plot(dt2_dt1p,rs1,'r');
hold on
plot(dt2_dt1,rs1,'b');
dips("Done");
