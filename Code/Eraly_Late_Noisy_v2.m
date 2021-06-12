%%% 
%%%for this file which is version 2 we listen to Dr Berber and add the
%%%interpolated noise to the signal.


tic;


clear;
clc;
%close all;

L = 16; 
sequence_length = L   ; %chaotic sequence Length;

SNR_dB = 8;           % Chip Eb/No


Rb = 200;  % bit rate = bits/sec
Rc = 100;   % chip rate = chips/bit   in other words spreading factor
Rs = 10;      % sample rate = samples/chip  

Rt = Rb*Rc*Rs;  %total rate



%What sort of sequence do you want correlated? Equal that to 1;
PN = 1;
CHAOS = 0;
NOISE = 0;
INTERP = 0;
ZOH = 1;

 x = zeros(1,L);                  % make the original x vector
 initial_condition = 0.1133;
 x(1) = initial_condition;
 for ii = 2:length(x),
    x(ii) = 1- 2*x(ii-1).^2;
 end;
 x_original = x;  % want to use this for correct noise calculation

%     %FOR CHAOS
%     x= zeros(1,L);                  % make the original x vector
%      initial_condition = 0.1133;
%      x(1) = initial_condition;
%      for ii = 2:length(x),
%         x(ii) = 1- 2*x(ii-1).^2;
%      end;


%%%================================NOISE SEQUENCES=========================
if PN ==1
    %FOR PN
    LFSR_connections = 7;
    sequence_type = 1;
    
    %x1 is the original PN vector
    [x1] = m_sequence_generator(LFSR_connections, sequence_type);

    x_original = x1;  % want to use this for correct noise calculation
    
    if ZOH ==1
       
        x = zeros(1,length(x1)*Rs); % Introduce ZOH sampling of the signal 
        jj = 1;
        for mm = 1:length(x1)
           for nn = 1:Rs
               x(jj) = x1(mm);
               jj = jj+1;
           end
            nn = 1; 
        end;
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
        x = y_zoh;
        
%         for(xx = 1:length(x)),        %%for coding purposes
%           
%             if x(xx)>=0
%                 x(xx) = 1;
%             else
%                 x(xx) = -1;
%             end
%         end
        
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
   noisy_x = noiseless_x + noise;
   
   
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

if PN ==1
    R = R.*(Rs./((2^LFSR_connections)-1))./Rc;
end

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

Rd = temp_d.^2;
Ra = temp_a.^2;
% Rd = temp_d;
% Ra = temp_a;

%%%========================================================================
%%%===========================error calculation============================
%%%========================================================================
error = Rd - Ra;





%%%========================================================================
%%%=================================GRAPHICS===============================
%%%========================================================================


%%% ============================ for plotting ZOH or interpolated signals
%%% and the original one against actual time. We superimpose them to make
%%% the visualization more meaningfull. 

time = 0:(1/Rt):(1/Rt)*(L*Rs-1);
x_original = upsample(x_original,Rs);
% figure1 = figure;

%  axes('Parent',figure1,'FontSize',18,'FontName','Times New Roman');
%  box('on');
%  hold('all');
% 
% 
% stem(time,x,'r.','DisplayName','Sampled Signal');
% hold on;
% stem(time,x_original,'MarkerSize',3,'MarkerFaceColor','blue','Marker','square','DisplayName','Original Signal');
% hold on;
% if NOISE ==1
%     stem(time,noisy_x,'g.','DisplayName','Noisy Vector');
% end
% if CHAOS == 1
% title(['Sampled & Original Chaotic Seq. with ' ,int2str(Rs),' Samples/Chip and ',int2str(L),' chips'],'fontsize',10);
% end
% 
% if PN ==1
% title(['Sampled and Original PN Seq. with ' ,int2str(Rs),' Samples/Chip and ',int2str(L),' chips'],'fontsize',10);
% end

% xlabel('Time (s)','FontSize',30,'FontName','Times New Roman','interpreter','latex');
% ylabel('Amplitude (V)','FontSize',30,'FontName','Times New Roman','interpreter','latex');
% 
% 
% legend('show');
% a = legend;
% set(a,'interpreter','latex');

%  axis([0 0.001  -1.2 1.2])
% xlabel({['Actual time based on Rb=',int2str(Rb),'bits/sec Rc=', int2str(Rc), 'and Rs']},'fontsize',12)
% ylabel('Amplitude','fontsize',12);
% ==========================================================================

% 
% 
%   % Create  ACF plot 
figure2 = figure;
if PN ==1
 axes('Parent',figure2,'FontSize',18,'FontName','Times New Roman');
box('on');
hold('all');
 

%plot(lags./(Rs),R(:,2)./max(R(:,2)),'k','DisplayName','On time','linewidth',2);
 hold on;

% Create plot
 plot(lags./(Rs), Rd, 'ko','DisplayName','Delayed','MarkerSize',5,'LineStyle', '-','linewidth',2);

hold on;
% Create plot
plot(lags./(Rs), Ra,'kx','MarkerSize',7,'DisplayName','Advanced','LineStyle', '-','linewidth',2);

hold on;
plot(lags./(Rs), error,'k','DisplayName','Error Signal','LineStyle', '--','linewidth',2);
end

if CHAOS ==1
    
     axes('Parent',figure2,'FontSize',18,'FontName','Times New Roman');
box('on');
hold('all');

plot(lags./(Rs), Rd./max(Rd), 'ko','DisplayName','Delayed','MarkerSize',5,'LineStyle', '-','linewidth',2);
plot(lags./(Rs), Ra./max(Ra),'kx','MarkerSize',7,'DisplayName','Advanced','LineStyle', '-','linewidth',2);
plot(lags./(Rs), error./max(error),'k','DisplayName','Error Signal','LineStyle', '--','linewidth',2);


    
end
% Create legend

%legend('FontSize',16);


xlabel('$$\delta$$ (in $$T_c$$)','fontsize',30,'FontName','Times New Roman','interpreter','latex');
ylabel('Normalized D($$\delta$$)','FontSize',30,'FontName','Times New Roman','interpreter','latex');
legend('show');

a = legend;
set(a,'interpreter','latex');
 
% % plot(lags,R(:,2)./max(R(:,2)))
% % hold on;
% % plot(lags,error./max(error),'r')
% % hold on;
% 
%  axis([-2.5 2.5 -1.5 1.5])

%%------------------------------------------------------------------

% if CHAOS == 1
% 
%     title(['R(\tau) , Chaotic Early, Late and on time with ' ,int2str(Rs),' Samples/Chip ',int2str(L),'chips',' SNR=',int2str(SNR_dB),'dB'],'fontsize',16);
%     %title(['Note : ',int2str(nb_points),' points from ',int2str(nb_gaussian),' gaussian samples   ;  Freq : F=',int2str(F),' Hz ; noise_power : \sigma^2 = ',num2str(noise_power),' W'],'FontSize',10);
% end;
% 
% if PN == 1
%     title(['R(\tau) , PN Early, Late and on time with ' ,int2str(Rs),' Samples/Chip ',int2str(L),'chips',' SNR=',int2str(SNR_dB),'dB'],'fontsize',16);
% end;




 toc;
%Elapsed_Time= t/60

