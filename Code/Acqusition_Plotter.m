%-----------------------------------------------------------------------%
%   NAME:   Acqusition Plotter                                           %
%   TITLE:  Plots acqusition ( cyclic searching) for the report purposes      %
%   AUTHOR: Ramin Vali                                 % 
%                                                                       %
%   DESCRIPTION: It also is used as the results.        %
%                                                                       %
%   DATE: 19/3/2008                                                      %
%                                                                       %
%   LAST MODIFICATION: Creation (RV)                                    %
%                                                                       %
%   OPERATING CONDITIONS: Normal operating conditions apply             %
%                                                                       %
%-----------------------------------------------------------------------%

clear;
clc;
tic;

%======================= VARIABLES TO CHANGE ============================

Z_threshold = 0;        % Inital Z value
step = 10;              % Step-size: again a balance between speed and accuracy
period = 1000;          % Period of the chaotic pilot
corr_period = 200;      % Correlation period
initial_condition = 0.2351;
Offset = 63;            % Chip offset (delay between rx and tx)
users = 10;              % Additional users
SNR_dB = 8;           % Chip Eb/No

Rs = 1;                 %Number of Samples
PN = 0;
CHAOS = 1;
ZOH =0;
INTERP =1;
%========================================================================


Offset = Offset*Rs;
% Set up the tx and rx arrays
rx_ptr = 0;             
tx_ptr = Offset; 


    tx_array = zeros(1, period);
    tx_array(1,1) = initial_condition;

if CHAOS == 1,
    for ii = 1:(period-1),
        tx_array(1,(ii+1)) = 1 - 2*tx_array(1,ii).^2;
    end;
    
elseif PN == 1,
    
    LFSR_connections = input ('please eneter the number if LFSR taps (2^x)-1 >= 10  ');
    sequence_type = 1;
    
    %x1 is the original PN vector
    [tx_array] = m_sequence_generator(LFSR_connections, sequence_type);

    x_original = tx_array;  % want to use this for correct noise calculation
        
end

if ZOH ==1,
    sequence = 0,
    tx_array = ZOHSample(tx_array,Rs,sequence);
    
elseif INTERP ==1,
    tx_array = interp(tx_array,Rs);
end;



rx_array = tx_array;

% Set up extra users
% if (users == 0),
%     IUI = 0;
% else
%     for ii = 1:users,
%         IUI(ii,1) = rand;
%     end;
% end;

IUI = zeros(users,corr_period);

%------------------------------ NOISE CALCULATIONS -----------------------
fs = 8000;          % sampling frequency (cancels)
sf = 100;           % Spreading factor

% Calculate the noise - do the long way
Pavg = sum(tx_array.^2) / (length(tx_array));   % Variance of the logistic map = 0.5
% Eb is the total energy of the signal divided by the number of bits:
Etot = Pavg .* (length(tx_array)/fs);
Eb = Etot./(length(tx_array)./sf);
% and Ec is the total energy of the signal divided by the number of chips
Ec = Etot./(length(tx_array));
SNR_lin = 10.^(SNR_dB/10);
No = Eb/SNR_lin;
sigma = No*fs/2;

%------------------ Find the Auto-correlation of the pilot ----------------

for jj = 1:period*Rs,
    % Take the next T samples:
    for ii = 1:corr_period*Rs-1,
        % Modulo the long way since matlab has no concept of A(0)
        if rx_ptr >= period*Rs,
            rx_ptr = 0;
        end;
        if tx_ptr >= period*Rs,
            tx_ptr = 0;
        end;
        % incriment pointers after doing the modulo 
        rx_ptr = rx_ptr + 1 ;
        tx_ptr = tx_ptr + 1;

        % Take the correlation samples from the fixed arrays using the
        % offset pointers
        rx_corr(ii) = rx_array(rx_ptr);
        tx_corr(ii) = tx_array(tx_ptr);   % add noise to the tx
    end;
    % Find the correlation sum
    R(jj) = (sum(rx_corr .* tx_corr)); 
    Z(jj) = .8;
    if rx_ptr >= period*Rs,
         rx_ptr = 0;
    else 
         rx_ptr = rx_ptr + 1;
    end;
end;


 
% Having found the correlation plot add the noise part:
% Repeat this N times for each threshold:


upper_ptr = 200;
    lower_ptr = 1;

    
    
    for jj = 1:period,
        
        % Generate extra users if required
        if (users == 0),
            total_IUI = 0;
        else
            for kk = 1:users
                for ii = 1:(200),
                    % Users transmit from the logistic map
                    if(ii == 1),
                        random_seed = rand;
                        while (random_seed == 1),
                            random_seed = rand;
                        end;                        
                        IUI(kk,ii) = random_seed;
                    else                        
                        IUI(kk,ii) = 1-2*(IUI(kk,ii-1)).^2;
                        if (IUI(kk,ii) > 0.9999),
                            while (IUI(kk,ii) > 0.9999),
                                IUI(kk,ii) = rand;
                            end; 
                        end;
                    end;
                end;
            end;
            
            % Total IUI is the sum of all users
            if (users == 1),
                total_IUI = IUI;
            else
                total_IUI = sum(IUI);
            end;
        end;
        
        % From mathematics:
        % Seperate the summation of the product into 2 components:
        %   1. The auto-correlation (already done)
        %   2. local-sequences*(Summation of users(done) + noise)
        % The final Z is the summation of both
       % noisy_tx = tx_array(lower_ptr:upper_ptr); %+ (randn(1,period/5)*sqrt(sigma) + total_IUI);
        
        % Move pointers
        upper_ptr = upper_ptr + 200;
        lower_ptr = lower_ptr + 200;
        if upper_ptr > 1000,
            upper_ptr = 200;
            lower_ptr = 1;
        end;
    
      
    end;
    figure1 = figure;
    axes('Parent',figure1,'FontSize',16,'FontName','Times New Roman');
box('on');
hold('all');
    
   plot(R./max(R));
   hold on;
   plot(Z, 'r', 'LineWidth',3)
   

   
  %axes('Parent',current,'FontSize',16);

   
   % Create textbox
annotation(figure1,'textbox','interpreter','latex','String',{'Threshold'},...
    'HorizontalAlignment','center',...
    'FontSize',22,...
    'FontName','Times New Roman',...
    'FitHeightToText','off',...
    'LineStyle','none',...
    'Position',[0.7266 0.7600 0.1656 0.08214]);
    
xlabel('Propagation delay in the channel (in $$T_c$$)','FontSize',30,'FontName','Times New Roman','interpreter','latex');
ylabel('Normalised $$Z$$ value','FontSize',30,'FontName','Times New Roman','interpreter','latex');

% Graphics

    
toc; 