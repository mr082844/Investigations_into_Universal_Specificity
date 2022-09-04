%% can traveling twin ever become older?
%close all
% functions
gamma = @(v) 1./sqrt(1-v.^2);
UIF_vel = @(v1_in,v2_in) (v1_in+v2_in)/(1 + v1_in*v2_in);
rel_vel = @(vUIF_in,v2_in) (v2_in-vUIF_in)/(v2_in*vUIF_in-1);

%% twins start from shared accelerated Frame

% given (all in UIF's perspective)
v1_1 = .4;   % 1st twin's travel speed as v/c (as seen by UIF)
v1_2 = .4;   % 1st twin's travel speed as v/c (as seen by UIF)
v2_1 = 0;  % second twin's travel speed period dt1 as v/c (as seen by UIF)
v2_2 = UIF_vel(v1_1,v1_1); % second twin's travel speed period dt2 as v/c (as seen by UIF)
dtp_1 = 5; % seconds after 1st twin starts travel until dt2 starts starts (as seen by UIF)

% determined (all in UIF's perspective)
x1_1 = v1_1 * dtp_1;       % distance 1st twin travels from t0' to t1' (as seen by UIF)
x2_1 = v2_1 * dtp_1;       % distance 2nd twin travels from t0' to t1' (as seen by UIF)
dx_1 = x1_1 - x2_1;        % distance between twins at t1' (as seen by UIF)
dv12_1_to_2 = v2_2 - v1_2; % relative velocity between twins from t1' to t2' (as seen by UIF)
dtp_2 = dx_1 / dv12_1_to_2;% time it takes for twins to rendezvous (as seen by UIF)
total_dtp = dtp_1 + dtp_2; % total time elapsed for experiment (as seen by UIF)
x1_2 = v1_2 * dtp_2;       % distance 1st twin travels from t1' to t2' (as seen by UIF)
x2_2 = v2_2 * dtp_2;       % distance 2nd twin travels from t1' to t2' (as seen by UIF)
total_dxp = x1_1+x1_2;

% twins time differentials (as seen by UIF)
dt1_dtp_1 = 1 / gamma(v1_1); % 1st twin's dt/dt' at t1' and t2' (as seen by UIF)
dt1_dtp_2 = 1 / gamma(v1_2); % 1st twin's dt/dt' at t1' and t2' (as seen by UIF)
dt2_dtp_1 = 1 / gamma(v2_1); % 2nd twin's dt/dt' at t1' (as seen by UIF)
dt2_dtp_2 = 1 / gamma(v2_2); % 2nd twin's dt/dt' at t1' (as seen by UIF)

% twins aged by:
age1_1 = dt1_dtp_1 * dtp_1;   % 1st twin's change age at t1'
age1_2 = dt1_dtp_2 * dtp_2;   % 1st twin's change age at t2'
total_age1 = age1_1 + age1_2; % 2nd twin's total change age for entire trip
age2_1 = dt2_dtp_1 * dtp_1;   % 2nd twin's change age at t1'
age2_2 = dt2_dtp_2 * dtp_2;   % 2nd twin's change age at t2'
total_age2 = age2_1 + age2_2; % 2nd twin's total change age for entire trip

% determine light tx/receive to hit at turn around
mid_t = dtp_1;               % time of turn around point (as seen by UIF)
mid_x1 = x1_1;               % dist of 1st twin at t1' (as seen by UIF)
mid_x2 = x2_1;               % dist of 2nd twin at t1' (as seen by UIF)
slope1_1 = 1/v1_1;           % time-space slope 1st twin from t0' to t1' (as seen by UIF)
yint1_1 = mid_t - slope1_1 * mid_x1;
slope1_2 = 1/v1_2;           % time-space slope 1st twin from t1' to t2' (as seen by UIF)
yint1_2 = mid_t - slope1_2 * mid_x1;
slope_light_tx = -1;         % time-space slope of light emission (as seen by UIF)
yint_light_tx = mid_t - slope_light_tx * mid_x2;
slope_light_rx = 1;          % time-space slope of light return trip (as seen by UIF)
yint_light_rx = mid_t - slope_light_rx * mid_x2;
if isinf(slope1_1)
    UIF_t_tx = mid_t - mid_x2*slope_light_tx;
    UIF_x_tx = 0;
else
    M_mat = [slope_light_tx -1; slope1_1 -1];
    b_mat = [-yint_light_tx; -yint1_1];
    xy = M_mat \ b_mat;
    UIF_x_tx = xy(1);
    UIF_t_tx = xy(2);
end
if isinf(slope1_2)
    UIF_t_rx = mid_t + (total_dxp-mid_x2)*slope_light_rx;
    UIF_x_rx = total_dxp;
else
    M_mat = [slope_light_rx -1; slope1_2 -1];
    b_mat = [-yint_light_rx; -yint1_2];
    xy = M_mat \ b_mat;
    UIF_x_rx = xy(1);
    UIF_t_rx = xy(2);
end
t1_tx_a = UIF_t_tx * dt1_dtp_1;
t1_rx_a = UIF_t_rx * dt1_dtp_1;
slope_light_tx_check = (mid_t - UIF_t_tx) / (mid_x2 - UIF_x_tx);
slope_light_rx_check = (UIF_t_rx - mid_t) / (UIF_x_rx - mid_x2);

% plot results on space-time graph
fig = figure(1);
subplot(1,2,2);
hold off
plot([0 x1_1 total_dxp], [0 dtp_1 total_dtp],'r-o','LineWidth',2); % 1st twin's world line
hold on
plot([0 x2_1 total_dxp], [0 dtp_1 total_dtp],'b-o','LineWidth',2); % 2nd twin's world line
plot([UIF_x_tx mid_x2 UIF_x_rx],[UIF_t_tx mid_t UIF_t_rx],'k--o','LineWidth',2); % light path of laser to detect turn around point

% annotate age
h_arrow = annotation('textarrow','String',{'"Stationary" Twin''s','Elapsed Time:'...
    ,sprintf('%0.2f [Seconds]',age1_1)},'HorizontalAlignment','left','FontSize',14);
h_arrow.Parent = fig.CurrentAxes;
h_arrow.X = [x1_1+dtp_1/6 x1_1];
h_arrow.Y = [dtp_1-dtp_1/12 dtp_1];
h_arrow = annotation('textarrow','String',{'"Stationary" Twin''s','Elapsed Time:'...
    ,sprintf('%0.2f [Seconds]',total_age1)},'HorizontalAlignment','left','FontSize',14);
h_arrow.Parent = fig.CurrentAxes;
h_arrow.X = [total_dxp+dtp_1/6 total_dxp];
h_arrow.Y = [total_dtp-dtp_1/12 total_dtp];
h_arrow = annotation('textarrow','String',{'"Traveling" Twin''s','Elapsed Time:'...
    ,sprintf('%0.2f [Seconds]',age2_1)},'HorizontalAlignment','right','FontSize',14);
h_arrow.Parent = fig.CurrentAxes;
h_arrow.X = [x2_1-dtp_1/6 x2_1];
h_arrow.Y = [dtp_1-dtp_1/12 dtp_1];
h_arrow = annotation('textarrow','String',{'"Traveling" Twin''s','Elapsed Time:'...
    ,sprintf('%0.2f [Seconds]',total_age2)},'HorizontalAlignment','right','FontSize',14);
h_arrow.Parent = fig.CurrentAxes;
h_arrow.X = [total_dxp-dtp_1/6 total_dxp];
h_arrow.Y = [total_dtp-dtp_1/12 total_dtp];

% annotate speed and relative speed
vrel12_1 = rel_vel(v1_1,v2_1);
h_text = text(x1_1/2+.1,dtp_1/2,{sprintf('v_{UIF} = %0.2fc',v1_1),sprintf('v_{REL} = %0.2fc',vrel12_1)}...
    ,'HorizontalAlignment','left','VerticalAlignment','middle','FontSize',14);
h_text.Parent = fig.CurrentAxes;
h_text = text(x2_1/2-.1,dtp_1/2,{sprintf('v_{UIF} = %0.2fc',v2_1),sprintf('v_{REL} = %0.2fc',vrel12_1)}...
    ,'HorizontalAlignment','right','VerticalAlignment','middle','FontSize',14);
h_text.Parent = fig.CurrentAxes;
vrel12_2 = rel_vel(v1_2,v2_2);
h_text = text(x1_1 + x1_2/2+.1,dtp_1+dtp_2/2,{sprintf('v_{UIF} = %0.2fc',v1_2),sprintf('v_{REL} = %0.2fc',vrel12_2)}...
    ,'HorizontalAlignment','left','VerticalAlignment','middle','FontSize',14);
h_text.Parent = fig.CurrentAxes;
h_text = text(x2_1+x2_2/2-.1,dtp_1+dtp_2/2,{sprintf('v_{UIF} = %0.2fc',v2_2),sprintf('v_{REL} = %0.2fc',vrel12_2)}...
    ,'HorizontalAlignment','right','VerticalAlignment','middle','FontSize',14);
h_text.Parent = fig.CurrentAxes;

% annotate light emit and receive times
h_arrow = annotation('textarrow','String',{'"Stationary" Twin''s',sprintf('Tx = %0.2f [Seconds]',t1_tx_a)},'HorizontalAlignment','left','FontSize',14);
h_arrow.Parent = fig.CurrentAxes;
h_arrow.X = [UIF_x_tx+3*dtp_1/6 UIF_x_tx];
h_arrow.Y = [UIF_t_tx-3*dtp_1/15 UIF_t_tx];
h_arrow = annotation('textarrow','String',{'"Stationary" Twin''s',sprintf('Rx = %0.2f [Seconds]',t1_rx_a)},'HorizontalAlignment','left','FontSize',14);
h_arrow.Parent = fig.CurrentAxes;
h_arrow.X = [UIF_x_rx+dtp_1/6 UIF_x_rx];
h_arrow.Y = [UIF_t_rx-dtp_1/6 UIF_t_rx];

% plot details
grid on
title({'Case 2: Twins Started Accelerated Relative to UIF'},'FontSize',16);
xlabel('Distance [Light-Second]','FontSize',16);
ylabel('Time [Second]','FontSize',16);
legend('"Sationary" Twin','"Traveling" Twin','Light Tx/Rx','Location','SE','FontSize',16);
axis equal

%% twins start from UIF

% given (all in UIF's perspective)
v1_1 = 0;   % 1st twin's travel speed as v/c (as seen by UIF)
v1_2 = 0;   % 1st twin's travel speed as v/c (as seen by UIF)
v2_1 = -vrel12_1;  % second twin's travel speed period dt1 as v/c (as seen by UIF)
v2_2 = -vrel12_2; % second twin's travel speed period dt2 as v/c (as seen by UIF)

% twins time differentials (as seen by UIF)
dt1_dtp_1 = 1 / gamma(v1_1); % 1st twin's dt/dt' at t1' and t2' (as seen by UIF)
dt1_dtp_2 = 1 / gamma(v1_2); % 1st twin's dt/dt' at t1' and t2' (as seen by UIF)
dt2_dtp_1 = 1 / gamma(v2_1);     % 2nd twin's dt/dt' at t1' (as seen by UIF)
dt2_dtp_2 = 1 / gamma(v2_2);     % 2nd twin's dt/dt' at t1' (as seen by UIF)

dtp_1 = age2_1/dt2_dtp_1; % seconds after 1st twin starts travel until dt2 starts starts (as seen by UIF)

% determined (all in UIF's perspective)
x1_1 = v1_1 * dtp_1;         % distance 1st twin travels from t0' to t1' (as seen by UIF)
x2_1 = v2_1 * dtp_1;       % distance 2nd twin travels from t0' to t1' (as seen by UIF)
dx_1 = x1_1 - x2_1;        % distance between twins at t1' (as seen by UIF)
dv12_1_to_2 = v2_2 - v1_2;   % relative velocity between twins from t1' to t2' (as seen by UIF)
dtp_2 = dx_1 / dv12_1_to_2;% time it takes for twins to rendezvous (as seen by UIF)
total_dtp = dtp_1 + dtp_2; % total time elapsed for experiment (as seen by UIF)
x1_2 = v1_2 * dtp_2;         % distance 1st twin travels from t1' to t2' (as seen by UIF)
x2_2 = v2_2 * dtp_2;       % distance 2nd twin travels from t1' to t2' (as seen by UIF)
total_dxp = x1_1+x1_2;

% twins aged by:
age1_1 = dt1_dtp_1 * dtp_1; % 1st twin's change age at t1'
age1_2 = dt1_dtp_2 * dtp_2; % 1st twin's change age at t2'
total_age1 = age1_1 + age1_2;     % 2nd twin's total change age for entire trip
age2_1 = dt2_dtp_1 * dtp_1;       % 2nd twin's change age at t1'
age2_2 = dt2_dtp_2 * dtp_2;       % 2nd twin's change age at t2'
total_age2 = age2_1 + age2_2;     % 2nd twin's total change age for entire trip

% determine light tx/receive to hit at turn around
mid_t = dtp_1;               % time of turn around point (as seen by UIF)
mid_x1 = x1_1;               % dist of 1st twin at t1' (as seen by UIF)
mid_x2 = x2_1;               % dist of 2nd twin at t1' (as seen by UIF)
slope1_1 = 1/v1_1;           % time-space slope 1st twin from t0' to t1' (as seen by UIF)
yint1_1 = mid_t - slope1_1 * mid_x1;
slope1_2 = 1/v1_2;           % time-space slope 1st twin from t1' to t2' (as seen by UIF)
yint1_2 = mid_t - slope1_2 * mid_x1;
slope_light_tx = -1;         % time-space slope of light emission (as seen by UIF)
yint_light_tx = mid_t - slope_light_tx * mid_x2;
slope_light_rx = 1;          % time-space slope of light return trip (as seen by UIF)
yint_light_rx = mid_t - slope_light_rx * mid_x2;
if isinf(slope1_1)
    UIF_t_tx = mid_t - mid_x2*slope_light_tx;
    UIF_x_tx = 0;
else
    M_mat = [slope_light_tx -1; slope1_1 -1];
    b_mat = [-yint_light_tx; -yint1_1];
    xy = M_mat \ b_mat;
    UIF_x_tx = xy(1);
    UIF_t_tx = xy(2);
end
if isinf(slope1_2)
    UIF_t_rx = mid_t - mid_x2*slope_light_rx;
    UIF_x_rx = 0;
else
    M_mat = [slope_light_rx -1; slope1_2 -1];
    b_mat = [-yint_light_rx; -yint1_2];
    xy = M_mat \ b_mat;
    UIF_x_rx = xy(1);
    UIF_t_rx = xy(2);
end
t1_tx = UIF_t_tx * dt1_dtp_1;
t1_rx = UIF_t_rx * dt1_dtp_1;
slope_light_tx_check = (mid_t - UIF_t_tx) / (mid_x2 - UIF_x_tx);
slope_light_rx_check = (UIF_t_rx - mid_t) / (UIF_x_rx - mid_x2);

% plot results on space-time graph
subplot(1,2,1);
hold off
plot([0 x1_1 total_dxp], [0 dtp_1 total_dtp],'r-o','LineWidth',2); % 1st twin's world line
hold on
plot([0 x2_1 total_dxp], [0 dtp_1 total_dtp],'b-o','LineWidth',2); % 2nd twin's world line
plot([UIF_x_tx mid_x2 UIF_x_rx],[UIF_t_tx mid_t UIF_t_rx],'k--o','LineWidth',2); % light path of laser to detect turn around point


% annotate age
h_arrow = annotation('textarrow','String',{'"Stationary" Twin''s','Elapsed Time:'...
    ,sprintf('%0.2f [Seconds]',age1_1)},'HorizontalAlignment','left','FontSize',14);
h_arrow.Parent = fig.CurrentAxes;
h_arrow.X = [x1_1+dtp_1/6 x1_1];
h_arrow.Y = [dtp_1-dtp_1/12 dtp_1];
h_arrow = annotation('textarrow','String',{'"Stationary" Twin''s','Elapsed Time:'...
    ,sprintf('%0.2f [Seconds]',total_age1)},'HorizontalAlignment','left','FontSize',14);
h_arrow.Parent = fig.CurrentAxes;
h_arrow.X = [total_dxp+dtp_1/6 total_dxp];
h_arrow.Y = [total_dtp-dtp_1/12 total_dtp];
h_arrow = annotation('textarrow','String',{'"Traveling" Twin''s','Elapsed Time:'...
    ,sprintf('%0.2f [Seconds]',age2_1)},'HorizontalAlignment','right','FontSize',14);
h_arrow.Parent = fig.CurrentAxes;
h_arrow.X = [x2_1-dtp_1/6 x2_1];
h_arrow.Y = [dtp_1-dtp_1/12 dtp_1];
h_arrow = annotation('textarrow','String',{'"Traveling" Twin''s','Elapsed Time:'...
    ,sprintf('%0.2f [Seconds]',total_age2)},'HorizontalAlignment','right','FontSize',14);
h_arrow.Parent = fig.CurrentAxes;
h_arrow.X = [total_dxp-dtp_1/6 total_dxp];
h_arrow.Y = [total_dtp-dtp_1/12 total_dtp];

% annotate speed and relative speed
vrel12_1 = rel_vel(v1_1,v2_1);
h_text = text(x1_1/2+.1,dtp_1/2,{sprintf('v_{UIF} = %0.2fc',v1_1),sprintf('v_{REL} = %0.2fc',vrel12_1)}...
    ,'HorizontalAlignment','left','VerticalAlignment','middle','FontSize',14);
h_text.Parent = fig.CurrentAxes;
h_text = text(x2_1-.1,dtp_1/2,{sprintf('v_{UIF} = %0.2fc',v2_1),sprintf('v_{REL} = %0.2fc',vrel12_1)}...
    ,'HorizontalAlignment','right','VerticalAlignment','middle','FontSize',14);
h_text.Parent = fig.CurrentAxes;
vrel12_2 = rel_vel(v1_2,v2_2);
h_text = text(x1_1 + x1_2/2+.1,dtp_1+dtp_2/2,{sprintf('v_{UIF} = %0.2fc',v1_2)...
    ,sprintf('v_{REL} = %0.2fc',vrel12_2)}...
    ,'HorizontalAlignment','left','VerticalAlignment','middle','FontSize',14);
h_text.Parent = fig.CurrentAxes;
h_text = text(x2_1+x2_2/2-.1,dtp_1+dtp_2/2,{sprintf('v_{UIF} = %0.2fc',v2_2)...
    ,sprintf('v_{REL} = %0.2fc',vrel12_2)}...
    ,'HorizontalAlignment','right','VerticalAlignment','middle','FontSize',14);
h_text.Parent = fig.CurrentAxes;

% annotate light emit and receive times
h_arrow = annotation('textarrow','String',{'"Stationary" Twin''s',sprintf('Tx = %0.2f [Seconds]',t1_tx)},'HorizontalAlignment','left','FontSize',14);
h_arrow.Parent = fig.CurrentAxes;
h_arrow.X = [UIF_x_tx+2*dtp_1/3 UIF_x_tx];
h_arrow.Y = [UIF_t_tx UIF_t_tx];
h_arrow = annotation('textarrow','String',{'"Stationary" Twin''s',sprintf('Rx = %0.2f [Seconds]',t1_rx)},'HorizontalAlignment','left','FontSize',14);
h_arrow.Parent = fig.CurrentAxes;
h_arrow.X = [UIF_x_rx+dtp_1/6 UIF_x_rx];
h_arrow.Y = [UIF_t_rx-dtp_1/6 UIF_t_rx];

% plot details
axis equal
grid on
title({'Case 1: Twins Started In UIF'},'FontSize',16);
xlabel('Distance [Light-Second]','FontSize',16);
ylabel('Time [Second]','FontSize',16);
legend('"Sationary" Twin','"Traveling" Twin','Location','SE','FontSize',16);
subplot(1,2,1)
sub121_xlim = xlim;
sub121_ylim = ylim;
subplot(1,2,2)
sub122_xlim = xlim;
sub122_ylim = ylim;
best_xlim = [floor(min(sub122_xlim(1), sub121_xlim(1))) ceil(max(sub122_xlim(2), sub121_xlim(2)))];
best_ylim = [floor(min(sub122_ylim(1), sub121_ylim(1))) ceil(max(sub122_ylim(2), sub121_ylim(2)))];
ylim(best_ylim);
yticks([best_ylim(1):best_ylim(2)]);
xlim(best_xlim)
xticks([best_xlim(1):best_xlim(2)]);
subplot(1,2,1)
ylim(best_ylim)
yticks([best_ylim(1):best_ylim(2)]);
xlim(best_xlim)
xticks([best_xlim(1):best_xlim(2)]);

close(fig);