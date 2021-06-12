% code for the synchronisation engine for multiple samples per chip. The
% noise componenet has to be sorted out. Basically we need to determine if
% the way the noise component is calculated when there are multiple samples
% per chip present. This is important. The average power
figure;
clc;

clear;

initial_condition = 0.231;

SNR_dB = 8;
% 
samples = 10;
sequence_length = 10;
spreading_factor = 100;
% Pavg = 0.5;       % average power of a chaotic sequence.
% 
% bit_rate = 200; % 200 bits persecond;
% chip_rate = bit_rate * spreading_factor;
% sample_rate = chip_rate * samples;
% 
% fs = chip_rate * samples;
% 
% % Eb is the total energy of the signal divided by the number of bits:
% Etot = Pavg .* (sequence_length/fs);
% Eb = (Etot./sequence_length./spreading_factor);
% % and Ec is the total energy of the signal divided by the number of chips
% Ec = Etot./(sequence_length);
% SNR_lin = 10.^(SNR_dB/10);
% No = Eb/SNR_lin;
% sigma = No*fs/2;



%%%Tryadding the noise to the chaotic sequence first then samping it in the
%%%square wave fasion. See what that gives you.!!

initial = zeros(1,sequence_length-1);   % this one is just for comarison reasons only.

x = zeros(samples,sequence_length-1); % to get and exact sequence length
                                      % we want, matlab arrays start at 1
x(:,1)=initial_condition;   %initial condition
initial(1,1) = initial_condition;
for (ii = 1:sequence_length-1),
    x(:,ii+1) = 1- 2.* (x(:,ii).^2); %makes a matrix of dimentions(sample,sequencelength)  the coloums contain the same numbers
    
    initial(1,ii+1) = 1-2.*(initial(1,ii).^2);  % we make this array( equal to the unsampled value for the x chaotic sequence,
                                                %for the matter of calculating the error of the  decimated values and the original sequence

end;                                 
 
    

    y  = [x(:,1).'];        % get the first colum and make it into the first row
    
    for (jj  = 2:sequence_length),
        
        y = [y,x(:,jj).'];      %keep appending to this row to make the transmission vector.
        
    end;
    
    % sorting the noise --------------------------------

        SNR_dB = 0;
    noise_var = 0.5/(spreading_factor*samples)./(10^(SNR_dB/10));  % sorts out the noise variance
    
    for kk = 1:length(y),
       y_noisy(kk) = y(kk) + sqrt(noise_var)*randn;
    
    end;
    
   % stem(y_noisy,'r');
    hold on;
    stem(y);
    
    % just trying the decimator on the y we have :
    y_decimate = decimate (y,samples,'fir');
   
    axis = 1:samples:samples*(length(x));
    %stem(axis,y_decimate,'k');
    
   % stem(axis,initial,'r');
    
    error1 = (y_decimate-initial);
   % stem(error1,'r');
    
    y_downsample = downsample(y,samples);
    
    stem(axis, y_downsample ,'r');
     error2 = (y_downsample-initial);
    %stem(error2);
    
    
    
%     %% sorting the IUI
%     
% %     users = 1;
% % %ISI initial conditions
% % for kk = 1:users,
% %     chaotic_interference(kk) = randn;
% % end;
% % 
% %  % ---------------- ISI
% %         for kk = 1:users,
% %             chaotic_interference(kk) = 1 - 2*chaotic_interference(kk).^2;
% %         end;
% %         
% %         total_ISI = sum(chaotic_interference);
