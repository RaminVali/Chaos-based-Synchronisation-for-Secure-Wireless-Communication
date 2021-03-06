%=========Noise and Chaos Energy calculations by Ramin Vali Date: 4/10/07


%% Cell 1 energy of chotic signals when and when not interpolated or ZOH
%%sampled
clc
clear

%To calculate the energy of a signal we can square the signal and multiply
%it by the time for the signal (\delta t) we need to calculate this \delta
%t hence we need to make some assumptions about the Rates at which we are
%working. There are three rates that we need to know about:
%   Rb = bit rate  and in bits/sec
%   Rc = the spreading factor and in chips/bit
%   Rs = sample number and in samples/chip
%   Hence the total rate is : Rt = Rb*Rc*Rs = Samples/sec

%the time between the samples is \delta t = 1/Rt.

%%%=========initialize variables:


Fin_time = 1; %finish time in secodns.

%we assume the following:
 
 Rb = 200; %That the bit rate is 200 bits per second
 Rc = 100; %That the spreading factor is 100 chips per bit
 Rs = 1;   %That there are 1 samples per chip
 
 %Hence we build our time vector which starts from zero and accomodates all
 %the values we generate:
 
 Rt = Rb*Rc*Rs;  %total rate
 
 time = 0:(1/Rt):Fin_time;
 time = time(1:length(time)-1);
 L = length(time);     %The length of the original chaotic vector length ( beofre interpolation)
sequence_length = L;
%first we generate Chaos vector
x = zeros(1,L);                  % make the original x vector
 initial_condition = 0.1133;
 x(1) = initial_condition;
 for ii = 2:length(x),
    x(ii) = 1- 2*x(ii-1).^2;
 end;

% stem(time,x) %debugging purposes only
 
 
 
 
 %second we calculate the energy of this initial vector, which is the
 %summation of the squared values times their duration. *this comes from
 %the Gilley J Transrypt international Inc.
 
 E_tot_orig = sum(abs(x).^2)/Rt;
 
 %%======================interpolation==================================
 
 Rs = 10;
 Rt = Rb*Rc*Rs;  %total rate
 
 time = 0:(1/Rt):Fin_time;
 time = time(1:length(time)-1);
 
 x_interpolated = interp(x,Rs);  %% generate interpolated signal
 
         
 E_tot_interp = sum(abs(x_interpolated).^2)/Rt;
 
 %figure; stem(time, x_interpolated)
 
 %=========================For ZOH signal======================
 
 %%%_____________________The ZOH signal___________________________
        
   Rs = 10;
   Rt = Rb*Rc*Rs;  %total rate
 time = 0:(1/Rt):Fin_time;
 time = time(1:length(time)-1);
 
 
        initial_condition = 0.1133;
         %% for the zero order hold sampling (square sampling method)
        x_zoh = zeros(Rs,sequence_length-1); % to get and exact sequence length
                                              % we want, matlab arrays start at 1
        x_zoh(:,1)=initial_condition;   %initial condition
        
        for (ii = 1:sequence_length-1),
            x_zoh(:,ii+1) = 1- 2.* (x_zoh(:,ii).^2); %makes a matrix of dimentions(sample,sequencelength)  the coloums contain the same numbers
        end;                                 
         
        y_zoh  = [x_zoh(:,1).'];        % get the first colum and make it into the first row
           
        
        for (jj  = 2:sequence_length),      
           y_zoh = [y_zoh,x_zoh(:,jj).'];      %keep appending to this row to make the transmission vector.      
        end;
        x_zoh = y_zoh;
        
        E_tot_ZOH = sum(abs(x_zoh).^2)/Rt;
    %stem(time,x_zoh)  %%debugging only
        status = 'finished cell 1'
 
%% Cell 2 Noise Energy when interpolated
clc

%we use the gaussian noise generator from matlab and we generate the noise power for Rs = 1, then we interpolate it! we can look at the variance  
%------------------------------ NOISE CALCULATIONS -----------------------
%%%=========initialize variables:
SNR_dB = 0;



%we assume the following:
 
 Rb = 200; %That the bit rate is 200 bits per second
 Rc = 100; %That the spreading factor is 100 chips per bit
 Rs = 1;   %That there are 1 samples per chip
 
  Rt = Rb*Rc*Rs;  %total rate
 
 time = 0:(1/Rt):Fin_time;
 time = time(1:length(time)-1);
 
 
    fs = Rt;          % sampling frequency (cancels)
    sf = Rc;           % Spreading factor

    % Calculate the noise - do the long way
    Pavg = sum(x.^2)/(length(x));   % Variance of the logistic map = 0.5
    % Eb is the total energy of the signal divided by the number of chips:
    Etot = sum(x.^2)/fs;
    Eb = Etot/((length(x)./(Rc)));
   
    SNR_lin = 10.^(SNR_dB/10);
    No = Eb/SNR_lin;  %noise PSD, watts per hertz
    sigma_squared = No*fs/2

noise_vector = (randn(1,length(x)).*sqrt(sigma_squared));

E_noise = sum(noise_vector.^2)/fs;

%  original signal
figure;stem(time,noise_vector ,'r','MarkerSize',4,'DisplayName','Noise Vector');
hold on;
stem (time,x,'MarkerSize',4,'DisplayName','Signal Vector');
title([' Original Chaotic Seq. with  ' ,int2str(SNR_dB),' dB of Noise & 1 S/C'],'fontsize',10);
xlabel({['Actual time based on Rb=',int2str(Rb),'bits/sec Rc=', int2str(Rc), 'and Rs']},'fontsize',12)
ylabel('Amplitude','fontsize',12);
legend('show');
%The noise is  now interpolated
%upsample the time vector fot graphing purposes
Rs = 100;   %That there are 1 samples per chip
 
  Rt = Rb*Rc*Rs;  %total rate
 
 time = 0:(1/Rt):Fin_time;
 time = time(1:length(time)-1);
 

noise_vector_interp = interp(noise_vector,Rs);
E_noise_interp = sum(noise_vector_interp.^2)/Rt;

figure;stem(time,noise_vector_interp ,'r','MarkerSize',4,'DisplayName','Noise Vector');
hold on;
stem (time,x_interpolated,'MarkerSize',4,'DisplayName','Signal Vector');
title([' interpolated Chaotic Seq. with  ' ,int2str(SNR_dB),' dB of Noise & ',int2str(Rs),' Sa/C'],'fontsize',10);
xlabel({['Actual time based on Rb=',int2str(Rb),'bits/sec Rc=', int2str(Rc), 'and Rs']},'fontsize',12)
ylabel('Amplitude','fontsize',12);
legend('show');




 





