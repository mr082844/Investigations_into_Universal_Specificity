function acceleration
close all
clear all
clc

% constants/variables
G = 6.6743e-11;
MS = 1.9891e30;
M=1e0*MS;
c = 3e8;
eT = c^2/2;
rs=G*M/eT;
 
% functions
GM = G*M;
SPE = @(r) GM./r;
SKE = @(v) v.^2/2;
gammaPE = @(r) 1./sqrt(1-SPE(r)/eT);
vel_free = @(r) sqrt(2*GM./r);
vel_spec2rel = @(r,v) gammaPE(r).*v;
gammaKE = @(v) 1./sqrt(1-SKE(v)/eT);
gammaT = @(r,v) gammaPE(r).*gammaKE(vel_spec2rel(r,v));
vel_spec = @(r) sqrt(2*GM./r)./gammaPE(r);
gammaT2 = @(r,v) 1./sqrt(1-(SKE(v)+SPE(r))/eT);
SIE = @(r,v) eT./(gammaT(r,v).^2);
acc_spec =  @(r) gammaPE(r).*GM.*sqrt(GM./r)./(sqrt(2)*eT*r.^2) - (GM)./(gammaPE(r).*sqrt(2).*r.^2.*sqrt(GM./r));
c0_c = @(r) 1./gammaT(r);
c0P_c = @(r) 1./gammaPE(r);
fg = @(r) GM./r.^2;
acc_spec2 = @(r,v) fg(r).*(1./(gammaPE(r).*gammaPE(r))...
    + gammaPE(r).^2.*(1-1./gammaKE(v).^2));

% range
r=flip([[rs:1:rs+1e5] [rs+1e5+1:1e5:rs+1e11] [rs+1e11+1:1e11:rs+1e15]]);
t_prime = time(r);
dt = t_prime(2:end)-t_prime(1:end-1);
t = [0 cumsum(dt).*gammaPE(r(2:end))];
dt = t(2:end)-t(1:end-1);
vel = vel_spec(r);
dv  = vel(2:end)-vel(1:end-1);
dv_dt = dv./dt;
dr = abs(r(2:end)-r(1:end-1));
dr_dt = dr./dt;
 
% Newtonian free-fall dynamics:https://mathpages.com/rr/s4-03/4-03.htm
 
% energies vs r
velp = vel_free(r);
eI = SIE(r,vel);
eK = SKE(vel);
eKp = SKE(vel_spec2rel(r,vel));
eP = SPE(r);
% eP = eT - eI - eK;
eT_sum = eI+eK+eP;

% internal, kinetic, and potential energy during freefall
figure(1);
hold off
h(2) = plot(r/rs,eI/eT,'LineWidth',2);
hold on
h(3) = plot(r/rs,eK/eT,'LineWidth',2);
h(4) = plot(r/rs,eP/eT,'LineWidth',2);
% h(1) =plot(r/rs,c0P_c(r),'k','LineWidth',2);
h(1) =plot(r/rs,eT_sum/eT,'--k','LineWidth',3);
h(5) =plot(r/rs,eK.*gammaPE(r).^2/eT,':k','LineWidth',1);
legend(h,'$\frac{e_I+\Delta e_K+\Delta e_P}{e_T}$','$\frac{e_I}{e_T}$'...
    ,'$\frac{\Delta e_K}{e_T}$','$\frac{\Delta e_P}{e_T}$'...
    ,'$\frac{\Delta e_{K|P}}{e_T}$','interpreter','latex','fontsize',20);
set(gca,'xscale','log')
grid on
xlabel('Distance From Center Mass [Schwarzschild Radii]'...
    ,'interpreter','latex','fontsize',20);
ylabel('Fraction of $e_T$ [-]','interpreter','latex','fontsize',20);
xlim([1e0 1e3]);
ylim([0 1]);
yticks([0:0.25:1]);
uistack(h(1:4),'top')
uistack(h(5),'top')
 
% acceleration during free fall
acc = acc_spec2(r,vel)/fg(rs);
figure(2)
plot(r/rs,acc);
set(gca,'xscale','log')
% set(gca,'yscale','log')
xlim([1e0 1e3]);
grid on


% net force=0 at rs=2
eKp(eKp>0.5*eT) = 0.5*eT;
eK = eKp./gammaPE(r).^2;
eI = eT - eK - eP;
figure (3)
hold off
h(2) = plot(r/rs,eI/eT,'LineWidth',2);
hold on
h(3) = plot(r/rs,eK/eT,'LineWidth',2);
h(4) = plot(r/rs,eP/eT,'LineWidth',2);
% h(1) =plot(r/rs,c0P_c(r),'k','LineWidth',2);
h(1) =plot(r/rs,eT_sum/eT,'--k','LineWidth',3);
h(5) =plot(r/rs,eK.*gammaPE(r).^2/eT,':k','LineWidth',1);
legend(h,'$\frac{e_I+\Delta e_K+\Delta e_P}{e_T}$','$\frac{e_I}{e_T}$'...
    ,'$\frac{\Delta e_K}{e_T}$','$\frac{\Delta e_P}{e_T}$'...
    ,'$\frac{\Delta e_{K|P}}{e_T}$','interpreter','latex','fontsize',20);
set(gca,'xscale','log')
grid on
xlabel('Distance From Center Mass [Schwarzschild Radii]'...
    ,'interpreter','latex','fontsize',20);
ylabel('Fraction of $e_T$ [-]','interpreter','latex','fontsize',20);
xlim([1e0 1e3]);
ylim([0 1]);
yticks([0:0.25:1]);
uistack(h(1:4),'top')
uistack(h(5),'top')

    function t = time(r)
        [r_sort,isort] = sort(r,'descend');
        [~,j_unsort] = sort(isort);
        a_sort = GM./r_sort.^2;
        a_sort(2:end) = geomean([a_sort(1:end-1);a_sort(2:end)],1);
        s_sort = max(r_sort) - r_sort;
        v0 = 0;
        t_sort = zeros(size(r));
        for ir = 2 : length(r)
            % 0 = at^2/2 + v0t - s
            this_a = -a_sort(ir)/2;
            this_b = v0;
            this_c = -(r_sort(ir) - r_sort(ir-1));
            this_t = (-this_b - sqrt(this_b^2-4*this_a*this_c))/(2*this_a);
            t_sort(ir) = this_t;
            v0 = a_sort(ir)*this_t;
        end
        % unsort
        t = t_sort(j_unsort);
    end
    function gm = geomean(X,dim)
        gm = (prod(X,dim)).^(1/size(X,dim));
    end
% nuclear forces
r_nuc = 0.8e-15;
m_nut = 1.674e-27;
rs_nut = G*m_nut/eT;
 
end
