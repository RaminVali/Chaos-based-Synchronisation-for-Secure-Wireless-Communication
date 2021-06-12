%%% Plots the SSB and DSB FFTs os ZOH and interpolated chaotic sequences.
%%%the SSB fft is done the matlab way. the DSB is done my own way.


tic;

clear;
clc;
%close all;

L = 1000; 
sequence_length = L   ; %chaotic sequence Length;
r = 10;
samples = r;



Rb = 100;  % bit rate = bits/sec
Rc = 100;   % chip rate = chips/bit   in other words spreading factor
Rs = 1;       % sample rate = samples/chip  

Rt = Rb*Rc*Rs;  %total rate

time = 0:(1/Rt):(1/Rt)*(L*Rs-1);  % time vector  Rs is 1 in here


x= zeros(1,L);                  % make the original x vector
initial_condition = 0.1133;
x(1) = initial_condition;
for ii = 2:length(x),
    x(ii) = 1- 2*x(ii-1).^2;
end;

%%%These are the calculations for single side band and double sideband 
%%%=======================================================================
%%%=================================DSB======================================
x_unsampled = x;
figure; stem(time,x_unsampled,'.');

title('Logistic Map With One Sample/Chip','fontsize',16)
xlabel('Time (s)')
ylabel('x(t)_unsampled')

y_unsampled = fft(x_unsampled);
y_unsampled = fftshift(y_unsampled);       %%shift the frequency of 0 to the center

f_Nyquist = Rt;  %Rt is the rate at which the signal is sampled so we assume it to be the nyquist rate.

F_max = f_Nyquist/2; %this is the maximum frquency resolvable

f_axis = linspace(-F_max, F_max, L);  
figure; plot(f_axis,((abs(y_unsampled)).^2)./F_max);   %same as multiplying something by its conjugate.

title('Double-Sided Power Spectrum of x(t)when unsampled')
xlabel('Frequency (Hz)')
ylabel('|X(f)|^2/Hz')

%%%===============================SSB=========================================

figure; %%_______ fft of the unsampled signal
NFFT = 2^nextpow2(L); % Next power of 2 from length of y
x_unsampled = x;
Y_unsampled = fft(x_unsampled,NFFT);
f_1 = (F_max)*linspace(0,1,NFFT/2);  %F_nyq is modified for Rs = 1;
%Plot single-sided amplitude spectrum.
plot(f_1,2*((abs(Y_unsampled(1:NFFT/2))).^2)/F_max)
%axis([0 50000 0 1.8]) %% arbitrary, just to sort out the numbers of the frequency axis
title('Single-Sided Power Spectrum of x(t)when unsampled')
xlabel('Frequency (Hz)')
ylabel('|X(f)|^2/Hz')



%%%========================================================================
%%%========================================================================
%%%_____________________The interpolated signal___________________________



%%%============================DSB=========================================
%%%===
%%% _________fft for the interpolated signal
x_interpolated = interp(x,r);  %% generate interpolated signal
  
Rs = r;       % sample rate = samples/chip  number of samples
Rt = Rb*Rc*Rs;  %total rate

time = 0:(1/Rt):(1/Rt)*(L*Rs-1);
figure;plot(time,x_interpolated,'.');
title(['Logistic Map With  in interpolation of ', int2str(r) ,' Samples/Chip']); %you need the [] brackets for the titlebar, otherwise it won't work.
xlabel('Time (s)')
ylabel('x(t)_interpolated')
  
y_interpolated = fft(x_interpolated);
y_interpolated = fftshift(y_interpolated);       %%shift the frequency of 0 to the center

f_Nyquist = Rt;  %Rt is the rate at which the signal is sampled so we assume it to be the nyquist rate.

F_max = f_Nyquist/2; %this is the maximum frquency resolvable

f_axis = linspace(-F_max, F_max, L*r);
figure; plot(f_axis,((abs(y_interpolated)).^2)./F_max);   %same as multiplying something by its conjugate.


title('Double-Sided Power Spectrum of x(t) when interpolatd')
xlabel('Frequency (Hz)')
ylabel('|X(f)|^2/Hz')

%%%================================SSB========================================
figure;%% _________fft for the interpolated signal
x_interpolated = interp(x,r);


NFFT = 2^nextpow2(L*r); %%adjusting the frequency component for the other one
f_2 = (F_max)*linspace(0,1,NFFT/2);


Y_interp = fft(x_interpolated,NFFT);

plot(f_2,2*(abs(Y_interp(1:NFFT/2)).^2)/F_max) %*2 because its the single sideband
title('Single-Sided Power Spectrum of x(t)when interpolatd')
xlabel('Frequency (Hz)')
ylabel('|X(f)|^2/Hz')

%%%========================================================================
%%%========================================================================
%%%========================================================================
%%%_____________________The ZOH signal___________________________


%%%===============================DSB=========================================
 %% for the zero order hold sampling (square sampling method)
x_zoh = zeros(samples,sequence_length-1); % to get and exact sequence length
                                      % we want, matlab arrays start at 1
x_zoh(:,1)=initial_condition;   %initial condition

for (ii = 1:sequence_length-1),
    x_zoh(:,ii+1) = 1- 2.* (x_zoh(:,ii).^2); %makes a matrix of dimentions(sample,sequencelength)  the coloums contain the same numbers
end;                                 
 
y_zoh  = [x_zoh(:,1).'];        % get the first colum and make it into the first row
   

for (jj  = 2:sequence_length),      
   y_zoh = [y_zoh,x_zoh(:,jj).'];      %keep appending to this row to make the transmission vector.      
end;
    % sample rate = samples/chip  number of samples
Rt = Rb*Rc*Rs;  %total rate

time = 0:(1/Rt):(1/Rt)*(L*Rs-1);
%figure;stem(time,y_zoh,'.');

% title(['Logistic Map With  ZOH sampling of ', int2str(r),' Samples/Chip and SF of ' , int2str(Rc)]);
% xlabel('Time (s)')
% ylabel('x(t)_ZOH')

  
Y_zoh = fft(y_zoh);
Y_zoh = fftshift(Y_zoh);       %%shift the frequency of 0 to the center

f_Nyquist = Rt;  %Rt is the rate at which the signal is sampled so we assume it to be the nyquist rate.

F_max = f_Nyquist/2; %this is the maximum frquency resolvable

f_axis = linspace(-F_max, F_max, L*r);
figure; plot(f_axis,((abs(Y_zoh)).^2)./F_max);   %same as multiplying something by its conjugate.we also divide by the nyquist frequency because its power spectral density.

title(['DSB-PSD of x(t) ZOH Sampled  ', int2str(r),' Samples/Chip and SF of ' , int2str(Rc)])
xlabel('Frequency (Hz)')
ylabel('|X(f)|^2/Hz')



%%%===============================SSB======================================
%%%===

figure; %% for the zero order hold sampling (square sampling method)

x_zoh = zeros(samples,sequence_length-1); % to get and exact sequence length
                                      % we want, matlab arrays start at 1
x_zoh(:,1)=initial_condition;   %initial condition
initial(1,1) = initial_condition;
for (ii = 1:sequence_length-1),
    x_zoh(:,ii+1) = 1- 2.* (x_zoh(:,ii).^2); %makes a matrix of dimentions(sample,sequencelength)  the coloums contain the same numbers

end;                                 
 
y_zoh  = [x_zoh(:,1).'];        % get the first colum and make it into the first row
   

for (jj  = 2:sequence_length),      
   y_zoh = [y_zoh,x_zoh(:,jj).'];      %keep appending to this row to make the transmission vector.      
end;
    
NFFT = 2^nextpow2(L*r); %%adjusting the frequency component for the other one
f_2 = (F_max)*linspace(0,1,NFFT/2);

Y_ZOH = fft(y_zoh,NFFT);

plot(f_2,2*(abs(Y_ZOH(1:NFFT/2)).^2)/F_max) 
title(['SSB-PSD of x(t) ZOH Sampled  ', int2str(r),' Samples/Chip and SF of ' , int2str(Rc)])
xlabel('Frequency (Hz)')
ylabel('|X(f)|^2/Hz')

toc;