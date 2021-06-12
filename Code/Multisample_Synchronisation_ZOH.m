%%% Multisample Synchronisation _ ZOH 

tic;

clear;
clc;
%close all;

initial_condition = 0.26354;

L = 1000; 
sequence_length = L   ; %chaotic sequence Length;
r = 10;
samples = r;



Rb = 100;  % bit rate = bits/sec
Rc = 100;   % chip rate = chips/bit   in other words spreading factor
Rs = 1;       % sample rate = samples/chip  

Rt = Rb*Rc*Rs;  %total rate

time = 0:(1/Rt):(1/Rt)*(L*Rs-1);  % time vector  Rs is 1 in here



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
  
Rs = r;       % sample rate = samples/chip  number of samples
Rt = Rb*Rc*Rs;  %total rate



%======================= VARIABLES TO CHANGE ============================
trials = 10^4;          % Balance between accuracy and speed: 10^3-10^5
Z_threshold = 0;        % Inital Z value
step = 10;              % Step-size: again a balance between speed and accuracy
period = sequence_length;          % Period of the chaotic pilot
corr_period = 200;      % Correlation period

Offset = 63*samples;            % Chip offset (delay between rx and tx)
%users = 8;              % Additional users
SNR_dB = 8;           % Chip Eb/No
%========================================================================

% Set up the tx and rx arrays
rx_ptr = 0;             
tx_ptr = Offset; 

tx_array = y_zoh;

rx_array = tx_array;
clear x_zoh;
clear y_zoh;


% Set up extra users-------------------------------------------------
% if (users == 0),
%     IUI = 0;
% else
%     for ii = 1:users,
%         IUI(ii,1) = rand;
%     end;
% end;

%IUI = zeros(users,corr_period);

%------------------------------ NOISE CALCULATIONS -----------------------
f_Nyquist = Rt;  %Rt is the rate at which the signal is sampled so we assume it to be the nyquist rate.
F_max = f_Nyquist/2; %this is the maximum frquency resolvable

fs = f_Nyquist;          % sampling frequency (cancels)
sf = Rc ;           % Spreading factor

% Calculate the noise - do the long way
Pavg = sum(tx_array.^2) / (length(tx_array));   % Variance of the logistic map = 0.5
% Eb is the total energy of the signal divided by the number of bits:
Etot = Pavg .* (length(tx_array)/fs);
Eb = Etot./(length(tx_array)./(sf*samples)); % because we need the number of bits which is the length of the sampled signal divided by the ( chips * samples) 
% and Ec is the total energy of the signal divided by the number of chips
Ec = Etot./(length(tx_array));
SNR_lin = 10.^(SNR_dB/10);
No = Eb/SNR_lin;
sigma_squared = No*fs/2;

% %------------------ Find the Auto-correlation of the pilot ----------------

for jj = 1:period*samples,
    % Take the next T samples:
    for ii = 1:corr_period*samples-1,
        % Modulo the long way since matlab has no concept of A(0)
        if rx_ptr >= period*samples,
            rx_ptr = 0;
        end;
        if tx_ptr >= period*samples,
            tx_ptr = 0;
        end;
        % incriment pointers after doing the modulo 
        rx_ptr = rx_ptr + 1;
        tx_ptr = tx_ptr + 1;

        % Take the correlation samples from the fixed arrays using the
        % offset pointers
        rx_corr(ii) = rx_array(rx_ptr);
        tx_corr(ii) = tx_array(tx_ptr);   % add noise to the tx
    end;
    % Find the correlation sum
    R(jj) = (sum(rx_corr .* tx_corr)); 

   
    
    if rx_ptr >= period*samples,
         rx_ptr = 0;
    else 
         rx_ptr = rx_ptr + 1;
    end;
end;

% %-------------------- Setup for the trials ------------------------------- 
% N_detect = 0;
% N_fail = 0;
% nn = 1;
% hh = waitbar(0, 'Estimated time remaining');
%  
% % Having found the correlation plot add the noise part:
% % Repeat this N times for each threshold:
% detections = zeros(1,100000);
% failures = zeros(1,100000);
% 
% for mm = 1:trials,
%     waitbar(mm/trials);
%     upper_ptr = 200;
%     lower_ptr = 1;
%     Z_threshold = 0;
%     
%     
%     for jj = 1:period,
%         
%         % Generate extra users if required
%         if (users == 0),
%             total_IUI = 0;
%         else
%             for kk = 1:users
%                 for ii = 1:(200),
%                     % Users transmit from the logistic map
%                     if(ii == 1),
%                         random_seed = 2*(rand-0.5);
%                         while (random_seed == 1),
%                             random_seed = 2*(rand-0.5);
%                         end;                        
%                         IUI(kk,ii) = random_seed;
%                     else                        
%                         IUI(kk,ii) = 1-2*(IUI(kk,ii-1)).^2;                        
%                         if (IUI(kk,ii) > 0.9999),
%                             while (IUI(kk,ii) > 0.9999),
%                                 IUI(kk,ii) = 2*(rand-0.5);
%                             end; 
%                         end;
%                     end;
%                     
%                     if IUI(kk,ii)>0,
%                       IUI_Binary(kk,ii)= 1;
%                     else
%                       IUI_Binary(kk,ii)= -1;  
%                     end;                   
%                 end;
%             end;
%             
%             % Total IUI is the sum of all users
%             if (users == 1),
%                 total_IUI = IUI_Binary;
%             else
%                 total_IUI = sum(IUI_Binary);
%             end;
%         end;
%         
%         % From mathematics:
%         % Seperate the summation of the product into 2 components:
%         %   1. The auto-correlation (already done)
%         %   2. local-sequences*(Summation of users(done) + noise)
%         % The final Z is the summation of both
%         noisy_tx = tx_array(lower_ptr:upper_ptr) .* (randn(1,period/5)*sqrt(sigma_squared) + total_IUI);
%         
%         % Move pointers
%         upper_ptr = upper_ptr + 200;
%         lower_ptr = lower_ptr + 200;
%         if upper_ptr > 1000,
%             upper_ptr = 200;
%             lower_ptr = 1;
%         end;
%     
%         N(jj) = sum(noisy_tx);  
%     end;
%     
%     % Compute Z
%     Z = (R + N).^2; 
%     
%     % Fast way to find the number of times over the threshold
%     index = find(Z > Z_threshold);
%     counter = 1;
%     
%     clear N_detect
%     clear N_fail;
%     
%     % Always find the max number of varaibles over the threshold
%     while(length(index) > 0),
%         index = find(Z > Z_threshold);
%         H1 = find(index == Offset+1);
%         N_detect(counter) = length(H1);
%         N_fail(counter) = length(index) - length(H1);
%         counter = counter + 1;
%         Z_threshold = Z_threshold + step;
%         Z_th_store(counter) = Z_threshold;
%     end;
%     
%     % Buffer up to 10000 since we can't add vectors of different lengths
%     N_detect(length(N_detect):100000) = 0;
%     N_fail(length(N_fail):100000) = 0;
%     detections = N_detect + detections;
%     failures = N_fail + failures;
%     n_store(mm,1:period) = N(1:period);
% end;
% 
% prob_failure = failures ./ max(failures);
% prob_detections = detections ./ max(detections);
% 
% prob_failure = prob_failure(1:length(Z_th_store));
% prob_detections = prob_detections(1:length(Z_th_store));
% 
% 
% 
% 
% MM = Z_th_store.';
% PP = prob_detections.';
% QQ = prob_failure.';
% 
% close(hh);
% 
% % % Graphics
% % plot(prob_failure, prob_detections,'b');
% % title('Receiver Operating Characteristic for Chaotic sequences');
% % xlabel('Probability of Failure');
% % ylabel('Probability of Detection');
% % grid on;
% % figure;
% % plot(Z_th_store, prob_detections(1:length(Z_th_store)), 'r');
% % hold on;
% % plot(Z_th_store, prob_failure(1:length(Z_th_store)), 'r');
% % grid on;
% % title('Cumulative Probability Density for Failure and Detection');
% % ylabel('Probability');
% % xlabel('Z-threshold');
%     
% toc; 