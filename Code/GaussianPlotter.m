%Simple Gussian plotter

sigma = 1;
mu = 0.0;
t = -6:0.001:6;
g = (1/(sqrt(2*pi)*sigma))*exp((-t.^2)/(2*sigma.^2));


plot(t,g./max(g));

