%-----------------------------------------------------------------------%
%   NAME:   CPSK                                                 	%
%   TITLE:  Simulation for the CPSK modulation scheme                   %
%   AUTHOR: Andrew Austin and Ramin Vali                                % 
%                                                                       %
%   DESCRIPTION: Thois Simulator simulates a CPSK system that uses a 
%   Sliding correlator rather than an integrating and dump circuit, The
%   gola is to add multiple samples per chip and find the correct way of
%   spreading the noise for the CPSK signal based on the BER verificiation
%   for multiple users.
%                                                                       %
%-----------------------------------------------------------------------%

clear;
clc;
tic;

for var = 1:3,

% --------------------- ONLY CHANGE THESE -----------------------------
% Produce a new BIT every tb seconds
fb = 100;          
N = fb;
% Produce SF chaotic samples per bit
SF = 100;             % Spreading factor

% Errors to be detected
err_detect = 10000;

if var ==1,
   Rs = 1
elseif var ==2,
    Rs = 5
   elseif var ==3,
    Rs = 10
end




users = 1;
% --------------------------------------------------------------------

% Initialize variables C style.
CPSK =  zeros (SF, fb);
CPSK_interp = zeros (SF*Rs, fb);
CPSK_trans =  zeros (SF, fb);
CPSK_trans_interp = zeros (SF*Rs, fb);
accumulator = zeros (1, fb);
input_vector = zeros (1, fb) ;
output_vector = zeros (1, fb);

noise =  zeros(SF*Rs, fb);


% Total number of sample points;
total_points = SF*fb;

% Rate of sampling
fs = SF*fb;          

SNR_dB = -2:1:8;                         % signal to noise ratio in decibels
SNR_lin = 10.^(0.1*SNR_dB);                 % signal to noise ratio in linear units

%Initialise variables
total_bit_errors = zeros(1,length(SNR_dB));
loop_counter = zeros(1,length(SNR_dB));

h = waitbar(0,'simulation running - please wait...');

for hh = 1:length(SNR_dB),
waitbar(hh/length(SNR_dB))    

while (total_bit_errors(hh) < err_detect),

%------------------ Produce Chaotic sequences ---------------------------
% Number of bits to be transmitted across
% Number of chaotic 'chips' per bit down 

for ii = 1:N,
    % Ensure the bit only changes every BIT period
    newBit = rand;
        
    if newBit < 0.5,
        input_vector(:,ii) = -1;
    else
        input_vector(:,ii) = 1;
    end;
        
    for kk = 1:SF,

        if (kk == 1),
            CPSK(kk,ii) = rand;
        elseif (kk <= SF),               
            CPSK(kk,ii) = 1 - 2*CPSK((kk-1) ,ii).^2;
        end;
            
        if newBit < 0.5,
            CPSK_trans(kk,ii) = CPSK(kk,ii) .* -1;
        else
            CPSK_trans(kk,ii) = CPSK(kk,ii);
        end;
                    
    end;      
    %%%Here we interpolate the transmitted signal
    
    CPSK_interp(:,ii) = interp(CPSK(:,ii),Rs);
    
    CPSK_trans_interp(:,ii) = interp(CPSK_trans(:,ii),Rs);
    
    
end;

%------------------------- IUI calculations ----------------------------
if users >1
    for uu = 1:users,
        for ii = 1:N,
            for kk = 1:SF,
                if (kk == 1),
                    CPSK(kk,ii) = rand;
                elseif (kk <= SF),               
                    CPSK_IUI(kk,ii) = 1 - 2*CPSK((kk-1) ,ii).^2;
                end;
            end;      
        end;
        CPSK_trans = CPSK_trans + CPSK_IUI;
    end;
else
    CPSK_trans = CPSK_trans;
end;
%------------------------- Noise calculations ----------------------------
% Calculate the energy per bit:
% Total energy / number of bits
% E_per_bit = sum(sum(CPSK.^2)) / (fb*fs);
% E_per_bit_store(hh) = E_per_bit;
P_avg = 0.5;  %for chaotic signals
E_per_bit = P_avg * (1/fb);

% Cross check - same formula only more complicated
%E_per_bit2 = sum(sum(CPSK.^2)) / (total_points * fb)
% Confirmed AA 26-03-06 14:06

% Find the 1-sided noise power spectral density
noise_PSD = E_per_bit / SNR_lin(hh);

% Noise BW is half the sampling rate
stand_dev = noise_PSD * fs / 2;



% randn = gaussian normal
noise = sqrt(stand_dev) * randn(SF,fb);


for ii = 1:N,
noise_interp(:,ii) = interp(noise(:,ii),Rs);
end;
% Add in the noise - matrix 
%CPSK_received = CPSK_trans.*fading + noise;  %nofading here!
%CPSK_received = CPSK_trans + noise; 
CPSK_received_interp = CPSK_trans_interp + noise_interp; 
clear noise_interp;
%-------------------------- END OF TRANSMITTER --------------------------

%-------------------------- START OF RECEIVER ---------------------------

for ii = 1:N,
    
    correlator = xcorr(CPSK_received_interp(:,ii),CPSK_interp(:,ii));

    
    if correlator(SF*Rs) > 0,
        output_vector(:,ii) = 1; 
    else
        output_vector(:,ii) = -1;
    end;
   
end;
%     CPSK_chaos_decoded = CPSK_received .* CPSK;
% 
%     accumulator = sum(CPSK_chaos_decoded);
%     
     for ii = 1:N, 
%         if accumulator(:,ii) < 0,
%             output_vector(:,ii) = -1;
%         else
%             output_vector(:,ii) = 1;
%         end;
        
        if output_vector(:,ii) ~= input_vector(:,ii),
            total_bit_errors(hh) = total_bit_errors(hh) + 1;
        else
            total_bit_errors(hh) = total_bit_errors(hh);
        end;
        
    end;        

% loop counter to keep track of how many times we have to run to get N
% errors
loop_counter(hh) = loop_counter(hh) + 1;

end; % end the while


end;
close(h) ;

%---------------------------- PLOTTING RESULTS ---------------------------

% BER is bits in error / total bits sent
BER = total_bit_errors ./ (fb*loop_counter);   

semilogy(SNR_dB, BER, 'ko-');
title('BER Plot for single user CPSK');
ylabel('Bit Error Rate');
xlabel('Eb/No (dB)');
hold on;
grid on;

for zz = 1:1:length(SNR_dB),
    
    psi = 0.5;      % variable used in BER eqn. its 0.5 all the time
    
    BER_eq(zz) = 0.5*erfc(( (psi/(SF/2) + ((users - 1)/(SF/2)) + ( SNR_lin(zz)).^-1)).^-.5);               %BER theory from CPSK 
    
end


semilogy(SNR_dB, BER_eq, 'rv-'); 
grid on;

legend('Simulated CPSK ','Theoretical CPSK - AWGN channel');
clear;
end;

toc;

