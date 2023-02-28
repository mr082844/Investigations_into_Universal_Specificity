c=1;
gamma = @(v) 1./sqrt(1-v.^2/c.^2); % v is fraction of c
avg_c = @(v) (c-v.*(v./c));
avg_c2 = @(v) c-v.*(2*c*v)./(c.^2+v.^2);
dv = 1e-4;
v=0:dv:(1-dv);
figure(1)
hold off
yyaxis left
plot(v,avg_c2(v),'r','linewidth',2)
ylim([0 1]);
yyaxis right
plot(v,avg_c2(v)-c./gamma(1-1./gamma(v)),'b','linewidth',2)
ylim([0 1]);
ax = gca;
ax.YAxis(1).Color ='k';
ax.YAxis(2).Color ='k';
grid on
% plot(v,(1-flip(avg_c(v,1).^2)),'r','linewidth',2)