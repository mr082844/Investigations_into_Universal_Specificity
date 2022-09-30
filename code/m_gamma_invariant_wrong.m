%% if KE derived from m_gamma being invariant (makes no sense at all)
c = 3e8;
gamma = @(v) 1 ./ sqrt(1 - v.^2/c^2);
KE = @(v) (1.5./gamma(v).^2 - .5).*.5.*v.^2; % resulting KE model if m_gamma is invariant
v = 0:1000:c;
figure(1);
plot(v/c,KE(v));