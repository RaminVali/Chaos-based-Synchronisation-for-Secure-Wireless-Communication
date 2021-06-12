%Third version, tidy up the code by using functions.
%put one chip and its corrresponding samples. introduce timing juitter in
%the channel and try to force the error to zero ( the REAL tracking loop!)
%Its still DDL though but good for a start


tic;


clear;
clc;
%close all;

L = 10; 
sequence_length = L   ; %chaotic sequence Length;

SNR_dB = 20;           % Chip Eb/No


Rb = 200;  % bit rate = bits/sec
Rc = 100;   % chip rate = chips/bit   in other words spreading factor
Rs = 1000;      % sample rate = samples/chip  

Rt = Rb*Rc*Rs;  %total rate



%What sort of sequence do you want correlated? Equal that to 1;
PN = 1;
CHAOS = 0;
NOISE =1;
INTERP = 0;
ZOH = 1;

 x = zeros(1,L);                  % make the original x vector
 initial_condition = 0.1133;
 x(1) = initial_condition;
 for ii = 2:length(x),
    x(ii) = 1- 2*x(ii-1).^2;
 end;
 x_original = x;  % want to use this for correct noise calculation


%%%================================PSEUDONOISE SEQUENCES=========================
if PN ==1
    %FOR PN
    LFSR_connections = 4;
    sequence_type = 1;
    
    %x1 is the original PN vector
    [x1] = m_sequence_generator(LFSR_connections, sequence_type);

    x_original = x1;  % want to use this for correct noise calculation
    
    if ZOH ==1
       
        sequence = 0;
        x = ZOHSample(x1,Rs,sequence); % Introduce ZOH sampling of the signal 0 for PN 1 for chaos
   
        L = 2^(LFSR_connections)-1;
        %% plotting purposes only
    %     Rt = Rb*Rc*Rs;  %total rate
    %    time = 0:(1/Rt):(1/Rt)*(L*Rs-1);
    %    figure; stem(time,x,'r.');
    end % end if ZOH ==1
    
    if INTERP ==1
        
        x_interpolated = interp(x1,Rs);  %% generate interpolated signal
         x = x_interpolated;
    end %end if INTERP ==1.
    
end; %end if  PN ==1

%%%================================CHAOTIC SEQUENCES=======================

if CHAOS ==1
    if ZOH == 1
        %%%========================================================================
        %%%========================================================================
        %%%========================================================================
        %%%_____________________The ZOH signal___________________________
        
  
       sequence = 1;
         x = ZOHSample(x1,Rs,sequence); % Introduce ZOH sampling of the signal 1 for chaos 0 for PN
        

        
            % sample rate = samples/chip  number of samples
        Rt = Rb*Rc*Rs;  %total rate
        
        time = 0:(1/Rt):(1/Rt)*(L*Rs-1);
        
         %figure; stem(time,x,'.');
        % 
        %  title(['Logistic Map With  ZOH sampling of ', int2str(r),' Samples/Chip and SF of ' , int2str(Rc)]);
        %  xlabel('Time (s)')
        %  ylabel('x(t)_ZOH')

        %=========================================================================
        %=======
        
    end  %end if ZOH = 1
    
   
    if INTERP ==1
         x = zeros(1,L);                  % make the original x vector
         initial_condition = 0.1133;
         x(1) = initial_condition;
         for ii = 2:length(x),
            x(ii) = 1- 2*x(ii-1).^2;
         end;
         x_original = x;  % want to use this for correct noise calculation
        x_interpolated = interp(x,Rs);  %% generate interpolated signal
         x = x_interpolated;
    end %end if INTERP ==1.
    
end %end if CHAOS ==1.  
    
    noiseless_x = x;
   
    
    
if NOISE == 1


    %------------------------------ NOISE CALCULATIONS -----------------------
    fs = Rt;          % sampling frequency (cancels)
    sf = Rc;           % Spreading factor

    % Calculate the noise - do the long way
    Pavg = sum(x.^2)/(length(x));   % Variance of the logistic map = 0.5
    % Eb is the total energy of the signal divided by the number of chips:
    Etot = sum(x.^2)/fs;
    Eb = Etot/((length(x)./(Rc)));
    
   
    % (length of x is the sample number so Es = 0.5 Ec)
   
    SNR_lin = 10.^(SNR_dB/10);
    No = Eb/SNR_lin;  %noise PSD, watts per hertz
    sigma_squared = No*fs/2
    
    if INTERP ==1
       noise =  interp((randn(1,length(x_original)).*sqrt(sigma_squared)),Rs);
    else
        noise =  (randn(1,length(x)).*sqrt(sigma_squared));
   end

    %-------------------------------------------------------------------------
  

   %initiate the noise vector
   noisy_x = noiseless_x; %+ noise;
   
   
  
   y_corr = [noisy_x ,noisy_x,noisy_x];  %%get ready for the correct correlation depending on prsence of noise
else
    sigma = 0;
    y_corr = [noiseless_x,noiseless_x,noiseless_x]; %%get ready for the correct correlation depending on prsence of noise
end; %endif NOISE =1;



%%% correlation

R = zeros(2*length(x)-1,3);


index = 1;
while(index <= 3)  %index is for the delays, There are 3 indecies, 1 for lag, 2 for on time and 3 for half a chip advance. 
    if index == 1
        y_ptr = 1;
    elseif index == 2
        y_ptr = (Rs/2)+ 1 ; 
    elseif index ==3
        y_ptr = Rs + 1 ;
    end
    ii = 1;
    while (y_ptr <(length(y_corr)-length(noiseless_x))),
        R(ii,index) = sum (noiseless_x(1:length(x)).*y_corr(y_ptr:y_ptr+length(x)-1));    
         y_ptr = y_ptr + 1; 
         ii = ii+1;
    end 
    index = index +1;
end
lags = (1:1:2*length(x)-1) - length(x); lags = (lags + ((Rs/2)-1)).'; %%from Salim's code. This is the x axis which is scaled for the different sampling sizes. 
% the goal is to make the  final x axis look like the x axis of one sample
% perchip. but just simply have mopre samples in the space between them.

% lags =  1:length(R(:,1));
% lags = lags.';


Rd = zeros(length(R(:,1)),1);
Ra = zeros(length(R(:,3)),1);

temp_d = R(:,1);  %R delay
temp_a = R(:,3);  %R advanced
% Rd(1) = temp_d(length(Rd));
% Ra(length(Ra)) = temp_a(1);
% for zz = 2:length(Rd)
%     Rd(zz) = temp_d(zz-1);
%     Ra(zz-1) = temp_a(zz);
%    Ra(zz) = temp_a(zz);
% end

Rd = temp_d;
Ra = temp_a;

%%%========================================================================
%%%===========================error calculation============================
%%%========================================================================
error = Rd - Ra;

error = error./max(error);

begin_error = (length(error)+1)/2-Rs+1;
end_error = (length(error)+1)/2 +1;

error = error.';

error_used = error(:,begin_error:end_error);


%%%========================================================================
%%%=================================GRAPHICS===============================
%%%========================================================================


%%% ============================ for plotting ZOH or interpolated signals
%%% and the original one against actual time. We superimpose them to make
%%% the visualization more meaningfull. 
% 
% time = 0:(1/Rt):(1/Rt)*(L*Rs-1);
% x_original = upsample(x_original,Rs);
% figure;stem(time,noisy_x,'g.','DisplayName','Noisy Vector');
%  hold on;
% stem(time,x,'r.','DisplayName','Sampled Vector');
% hold on;
% stem(time,x_original,'.','DisplayName','Original Vector');
% 
%  
% 
% if CHAOS == 1
% title(['Sampled & Original Chaotic Seq. with ' ,int2str(Rs),' Samples/Chip and ',int2str(L),' chips'],'fontsize',16);
% end
% 
% if PN ==1
% title(['Sampled and Original PN Seq. with ' ,int2str(Rs),' Samples/Chip and ',int2str(L),' chips'],'fontsize',16);
% end
% 
% legend('show');
% xlabel({['Actual time based on Rb=',int2str(Rb),'bits/sec Rc=', int2str(Rc), ', SNR = ',int2str(SNR_dB),'dB']},'fontsize',16)
% ylabel('Amplitude','fontsize',16);
% % ==========================================================================
% 
% % 
% % 
% %   % Create  ACF plot 
%   % Create plot
% figure; plot(lags./(Rs),R(:,2)./max(R(:,2)),'DisplayName','On time ACF','linewidth',1);
%  hold on;
% 
% % Create plot
% plot(lags./(Rs), Rd./max(Rd), 'DisplayName','Late','Color',[0 1 0],'linewidth',1);
% hold on;
% % Create plot
% plot(lags./(Rs), Ra./max(temp_a),'DisplayName','Early','Color',[1 0 0],'linewidth',1);
% hold on;
% plot(lags./(Rs), error./max(error),'DisplayName','e(\tau)','Color',[0 0 0],'linewidth',1);
% 
% 
% 
% % Create legend
% legend('show');
% xlabel('\tau (in T_c)','fontsize',16);
% ylabel('Normalized R(\tau)','fontsize',16);
% title(['R(\tau) , PN Early, Late and on time with ' ,int2str(Rs),' Samples/Chip ',int2str(L),'chips'],'fontsize',16);
% 
% % plot(lags,R(:,2))
% % hold on;
% % plot(lags,error,'r')
% % hold on;
% 
% 
% axis([-3 3 -1.5 1.5])
% 
% 
% if CHAOS == 1
% 
%     title(['R(\tau) , Chaotic Early, Late and on time with ' ,int2str(Rs),' Samples/Chip ',int2str(L),'chips',' SNR=',int2str(SNR_dB),'dB'],'fontsize',12);
%     %title(['Note : ',int2str(nb_points),' points from ',int2str(nb_gaussian),' gaussian samples   ;  Freq : F=',int2str(F),' Hz ; noise_power : \sigma^2 = ',num2str(noise_power),' W'],'FontSize',10);
% end;
% 
% if PN == 1
%     title(['R(\tau) , PN Early, Late and on time with ' ,int2str(Rs),' Samples/Chip ',int2str(L),'chips',' SNR=',int2str(SNR_dB),'dB'],'fontsize',12);
% end;




 toc;
%Elapsed_Time= t/60

