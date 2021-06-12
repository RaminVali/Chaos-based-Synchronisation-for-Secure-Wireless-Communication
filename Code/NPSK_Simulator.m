%-----------------------------------------------------------------------%
%   NAME:   NPSK                                                     	%
%   TITLE:  Simulation for the DCSK modulation scheme                   %
%   AUTHOR:  Ramin Vali                                                 % 
%                                                                       %
%   DESCRIPTION:       %
%                                                                       %
%                                                                       %
%   DATE:                                                     %
%                                                                       %
%   LAST MODIFICATION: 
%                                                                       %
%                                                                       %
%   OPERATING CONDITIONS: None                                          %
%                                                                       %
%-----------------------------------------------------------------------%

clear;
clc;
tic;
% --------------------- ONLY CHANGE THESE -----------------------------
% Produce a new BIT every tb seconds
fb = 10^2;          
N = fb;
% Produce SF chaotic samples per bit
SF = 100;             % Spreading factor

% Errors to be detected
err_detect = 10;

users = 1;
% --------------------------------------------------------------------

% Initialize variables C style.

NPSK = zeros (SF, fb);
accumulator = zeros (1, fb);
input_vector = zeros (1, fb) ;
output_vector = zeros (1, fb);

noise =  zeros(SF, fb);


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

        NPSK(kk,ii) = randn;
            
        if newBit < 0.5,
            NPSK_trans(kk,ii) = NPSK(kk,ii) .* -1;
        else
            NPSK_trans(kk,ii) = NPSK(kk,ii);
        end;
                    
    end;      
    
end;

%------------------------- IUI calculations ----------------------------
if users >1
    for uu = 1:users,
        for ii = 1:N,
            for kk = 1:SF,
                             
                    NPSK_IUI(kk,ii) = randn;
               
            end;      
        end;
        NPSK_trans = NPSK_trans + NPSK_IUI;
    end;
else
    NPSK_trans = NPSK_trans;
end;
%------------------------- Noise calculations ----------------------------
% Calculate the energy per bit:
% Total energy / number of bits
E_per_bit = sum(sum(NPSK.^2)) / (fb*fs);

% Cross check - same formula only more complicated
%E_per_bit2 = sum(sum(NPSK.^2)) / (total_points * fb)
% Confirmed AA 26-03-06 14:06

% Find the 1-sided noise power spectral density
noise_PSD = E_per_bit / SNR_lin(hh);

% Noise BW is half the sampling rate
stand_dev = noise_PSD * fs / 2;

% randn = gaussian normal
noise = sqrt(stand_dev) * randn(SF,fb);

% Add in the noise - matrix 
%NPSK_received = NPSK_trans.*fading + noise;  %nofading here!
NPSK_received = NPSK_trans + noise; 
%-------------------------- END OF TRANSMITTER --------------------------

%-------------------------- START OF RECEIVER ---------------------------


    NPSK_chaos_decoded = NPSK_received .* NPSK;

    accumulator = sum(NPSK_chaos_decoded);
    
    for ii = 1:N, 
        if accumulator(:,ii) < 0,
            output_vector(:,ii) = -1;
        else
            output_vector(:,ii) = 1;
        end;
        
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
title('BER Plot for single user NPSK');
ylabel('Bit Error Rate');
xlabel('Eb/No (dB)');
hold on;
grid on;

for zz = 1:1:length(SNR_dB),
    
    
    
    BER_eq(zz) = 0.5*erfc(((users+1)/(SF/2))+ ((SNR_lin(zz)).^-1).^-0.5);                                   %BER theory from NPSK 
    
end


semilogy(SNR_dB, BER_eq, 'rv-'); 
grid on;

legend('Simulated NPSK ','Theoretical NPSK - AWGN channel');

toc;

