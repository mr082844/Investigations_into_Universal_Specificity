c    = 299792458;  % [m/s] speed of light
v1   = 0.5*c;
v2   = 0;
x0_1 = 1*c;
x0_2 = 5*c;
t0   = 0;

% first pulse off of 2 from 1
dv_2   = c-v2;
dtp1_2 = (x0_2 - x0_1)/dv_2;
dxp1_2 = v2*dtp1_2;
xp1_2  = dxp1_2 + x0_2;
tp1_2  = t0 + dtp1_2;

% first pulse off of 1 from 2
dv_1   = -c-v1;
dtp1_1 = (x0_1 - x0_2)/dv_1;
dxp1_1 = v1*dtp1_1;
xp1_1  = dxp1_1 + x0_1;
tp1_1  = t0 + dtp1_1;

% first return pulse to 1
dt1_1 = 2*(x0_1 - x0_2)/dv_1;
dx1_1 = v1*dt1_1;
x1_1  = dx1_1 + x0_1;
t1    = t0 + dt1_1;

% first return pulse to 2
dt1_2 = 2*(x0_2 - xp1_1)/dv_2;
dx1_2 = v2*dt1_2;
x1_2  = dx1_2 + x0_2;
if 1e-12 < abs(dt1_2 - dt1_1)
    error("something is wrong");
end
% t1  = t0 + dt1_2; % dt1_2 should be the same as dt1_1

% second pulse off of 2 from 1
dtp2_2 = (x1_2 - x1_1)/dv_2;
dxp2_2 = v2*dtp2_2;
xp2_2  = dxp2_2 + x1_2;
tp2_2  = t1 + dtp2_2;

% second pulse off of 1 from 2
dtp2_1 = (x1_1 - x1_2)/dv_1;
dxp2_1 = v1*dtp2_1;
xp2_1  = dxp2_1 + x1_1;
tp2_1  = t1 + dtp2_1;

% second return pulse to 1
dt2_1 = 2*(x1_1 - x1_2)/dv_1;
dx2_1 = v1*dt2_1;
x2_1  = dx2_1 + x1_1;
t2    = t1 + dt2_1;

% second return pulse to 2
dt2_2 = 2*(x1_2 - xp2_1)/dv_2;
dx2_2 = v2*dt2_2;
x2_2  = dx2_2 + x1_2;
if 1e-12 < abs(dt2_2 - dt2_1)
    error("something is wrong");
end

% plot
hold off
plot([x0_1 xp1_1 x1_1 xp2_1 x2_1]/c,[t0 tp1_1 t1 tp2_1 t2],'r-x'...
    ,[x0_2 xp1_2 x1_2 xp2_2 x2_2]/c,[t0 tp1_2 t1 tp2_2 t2],'b-x','LineWidth',2);
xticks(0:(x0_2/c+1));
yticks(0:100);
grid on
xlabel('Space [Light-Seconds]')
ylabel('Time [Seconds]');
set(gca,'FontSize',14)
hold on
quiver(x0_2/c,t0,(xp1_1-x0_2)/c,tp1_1-t0,0,'--k');
quiver(x0_1/c,t0,(xp1_2-x0_1)/c,tp1_2-t0,0,'--k');
quiver(xp1_2/c,tp1_2,(x1_1-xp1_2)/c,t1-tp1_2,0,'--k');
quiver(xp1_1/c,tp1_1,(x1_2-xp1_1)/c,t1-tp1_1,0,'--k');
quiver(x1_2/c,t1,(xp2_1-x1_2)/c,tp2_1-t1,0,'--k');
quiver(x1_1/c,t1,(xp2_2-x1_1)/c,tp2_2-t1,0,'--k');
quiver(xp2_2/c,tp2_2,(x2_1-xp2_2)/c,t2-tp2_2,0,'--k');
quiver(xp2_1/c,tp2_1,(x2_2-xp2_1)/c,t2-tp2_1,0,'--k');
% axis equal
xlim([0 (x0_2/c+1)]);
% this_ylim =ylim;
ylim([0 floor(t2)+1])
