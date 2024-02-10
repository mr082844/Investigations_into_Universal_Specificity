G = 6.67430e-11;% [N m^2/kg^2] gravitational constant
c = 299792458;  % [m/s] speed of light
e_T = 0.5*c^2;  % [m^2/s^2] total specific energy
Ms = 1.9891e30; % [kg] mass of sun
MBH = 3.8*Ms;   % [kg] mass of smallest known black hole (2008): https://www.scientificamerican.com
Rs = G*MBH/e_T; % [m] radius of event horizon
gamma_g_SQ_inv = @(r,M) (1-(G*M./r)/e_T); % gravitational lorentz factor squared
x = [1:1e-2:1000];
r = x*Rs;
c_0_div_c_SQ = (gamma_g_SQ_inv(r,MBH));
e_P_div_e_T = (G*MBH./r)/e_T;
e_k_div_e_T_Free = e_P_div_e_T.*c_0_div_c_SQ;
figure(1)
hold off
plot(x,ones(size(x)),'k--','LineWidth',2);
hold on
plot(x,e_P_div_e_T,'r','LineWidth',2);
plot(x,e_k_div_e_T_Free,'g:','LineWidth',2);
plot(x,c_0_div_c_SQ,'b','LineWidth',2);
ylim([0 1]);
set(gca,'xscale','log')
legend({'$\frac{e_I}{e_T} + \frac{\Delta e_K}{e_T} + \frac{\Delta e_P}{e_T}$'...
    ,'$\frac{\Delta e_P}{e_T}=\frac{1}{R_s}$'...
    ,'$\frac{\Delta e_K}{e_T}$ in Free Fall'...
    ,'$\frac{e_I + \Delta e_K}{e_T}=\frac{c_I^2 + v^2}{c_\infty^2}=\frac{c_{R_s}^2}{c_\infty^2}$'}...
    ,'interpreter','latex'...
    ,'location','East');
xlabel('Radius to Mass [R_s]');
yticks([0:.1:1]);
grid on
ax = gca;
ax.FontSize = 16;