clear
inv_gamma_P = [0:1e-3:1];
inv_gamma_K = inv_gamma_P;
inv_gamma_P_mat = meshgrid(inv_gamma_P,inv_gamma_P);
inv_gamma_K_mat = meshgrid(inv_gamma_K,inv_gamma_K)';
inv_gamma_T_mat = (inv_gamma_P_mat.*inv_gamma_K_mat);
figure(1)
hold off
surf(inv_gamma_K_mat.^2,inv_gamma_P_mat.^2,inv_gamma_T_mat.^2,'EdgeColor','none');
hold on
plot3(inv_gamma_K,ones(size(inv_gamma_K)),inv_gamma_K,'k');
plot3(ones(size(inv_gamma_P)),(inv_gamma_P),(inv_gamma_P),'k');
plot3(zeros(size(inv_gamma_P)),(inv_gamma_P),zeros(size(inv_gamma_P)),'k');
plot3(inv_gamma_K.^2,zeros(size(inv_gamma_K)),zeros(size(inv_gamma_K)),'k');
xlabel('$\frac{1}{\gamma_K}$','Interpreter','latex','FontSize',20);
ylabel('$\frac{1}{\gamma_P}$','Interpreter','latex','FontSize',20);
zlabel('$\frac{1}{\gamma_T}$','Interpreter','latex','FontSize',20);
title('$\frac{1}{\gamma_T}=f(\frac{1}{\gamma_P},\frac{1}{\gamma_K})$','Interpreter','latex','FontSize',20);
grid on
colorbar
axis equal
