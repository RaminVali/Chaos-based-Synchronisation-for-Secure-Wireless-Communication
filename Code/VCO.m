%% VCO
clear;

pi = 3.1415;
Fc = 10; % in Hz
t = 0:0.001:1;

for ii = 1:length(t)
output(ii) = sin(2*pi*Fc*t(ii));
end;

plot(t,output);