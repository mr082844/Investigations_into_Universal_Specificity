c   = 299792458;  % [m/s] speed of light
G   = 6.6744e-11; % [m^3/(kg s)] gravitational constant
Me  = 5.97219e24; % [kg] earth's mass
M = Me;
e_T = 0.5*c^2;
r_s = G*M/e_T;
r_in_r_S = 1:1/100:15;
r = r_in_r_S*r_s;
c_0_c = sqrt(1-r_s./r);
e_I_e_T = c_0_c.^2;
del_e_I_e_T = (r_s)./(r).^2;
figure(1);
subplot(2,1,1);
hold off
plot(r_in_r_S,e_I_e_T,'b','LineWidth',2);
set(gca,'FontSize',14)
ylabel('$\frac{e_I}{e_T}$','Interpreter','latex','FontSize',16)
xticks(0:15);
yticks(0:.1:1);
grid on
subplot(2,1,2);
hold off
plot(r_in_r_S,del_e_I_e_T*r_s,'b','LineWidth',2);
set(gca,'FontSize',14)
ylabel('$\nabla\frac{e_I}{e_T}$','Interpreter','latex','FontSize',16)
xlabel('r_s','FontSize',16)
xticks(0:15);
yticks(0:.1:1);
grid on