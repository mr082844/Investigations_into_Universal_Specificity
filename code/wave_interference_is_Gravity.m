clear
close all
clc
t = 0:1e-8:1e-6;
lambd=1; % [hz]
c=3e6;
wave1 = @(t) cos(2*pi/lambd*(c*t));
wave2 = @(t) cos(pi/2 - 2*pi/lambd/sqrt(2)*(c*t));
figure(1)
hold off
plot(t,wave1(t))
hold on
plot(t,wave2(t))
plot(t,wave1(t)+wave2(t))
legend('wave1','wave2','wave1+wave2');