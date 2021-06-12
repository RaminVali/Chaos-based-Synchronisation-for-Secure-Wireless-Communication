%-----------------------------------------------------------------------%
%   NAME:   CPSK                                                 	%
%   TITLE:  Simulation for the DCSK modulation scheme                   %
%   AUTHOR: Andrew Austin and Ramin Vali                                % 
%                                                                       %
%   DESCRIPTION: Simulator for the DCSK modulation scheme in the        %
%           presence of AWGN and channel Fading (Rayleigh)              %
%           Produces neater curves as we run until x number of errors   %
%           have been detected                                          %
%                                                                       %
%           The results from running this simulator are shown in        %
%                                                                       %
%                                                                       %
%   DATE: 30-03-06                                                      %
%                                                                       %
%   LAST MODIFICATION: New modification to make the simulation time     %
%           shorter by initializing the variables (RV)                  % 
%           New modification to run until at least N errors have        %   
%           been detected, in this case fb is really irrelevent (AA, RV)%
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
fb = 10^3;          
N = fb;
% Produce SF chaotic samples per bit
SF = 100;             % Spreading factor

% Errors to be detected
err_detect = 100;

users = 1;
% --------------------------------------------------------------------

% Initialize variables C style.

CPSK = zeros (SF, 2*fb);
accumulator = zeros (1, 2*fb);
input_vector_uncoded = zeros(1,fb);
input_vector = zeros (1, 2*fb) ;
output_vector = zeros (1, 2*fb);

code = zeros(1,2*fb); %for coding purposes

noise =  zeros(SF, fb);


% Total number of sample points;
total_points = SF*fb;

% Rate of sampling
fs = SF*fb;          

SNR_dB = 2:1:7;                         % signal to noise ratio in decibels
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
        input_vector_uncoded(:,ii) = 0;
    else
        input_vector_uncoded(:,ii) = 1;
    end;
    
end;
    
%%%Jimmy's algorithm for coding=========================================
K=9;
trel=poly2trellis(K, [561 753]);
    
    
code=convenc(input_vector_uncoded,trel);  %input vector binary 0 and 1
%%%==================================================================
for jj = 1:2*N

    if code(:,jj) == 0,
            input_vector(:,jj) = -1;
        else
            input_vector(:,jj) = 1;
    end;

    
        
        for kk = 1:SF,

            if (kk == 1),
                CPSK(kk,jj) = rand;
            elseif (kk <= SF),               
                CPSK(kk,jj) = 1 - 2*CPSK((kk-1) ,jj).^2;
            end;

            if input_vector(:,jj) == -1,
                CPSK_trans(kk,jj) = CPSK(kk,jj) .* -1;
            else
                CPSK_trans(kk,jj) = CPSK(kk,jj);
            end;

        end;      
    
   
end;
%------------------------- IUI calculations ----------------------------
if users >1
    for uu = 1:users,
        for ii = 1:2*N,
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
E_per_symbol = sum(sum(CPSK.^2)) / (fb*fs);

% Cross check - same formula only more complicated
%E_per_bit2 = sum(sum(CPSK.^2)) / (total_points * fb)
% Confirmed AA 26-03-06 14:06

%Eb = 2*Es
% Find the 1-sided noise power spectral density
noise_PSD = 2*E_per_symbol / SNR_lin(hh);  %For the symbol

% Noise BW is half the sampling rate
stand_dev = noise_PSD * fs / 2;

% randn = gaussian normal
noise = sqrt(stand_dev).* randn(SF,2*fb);

% Add in the noise - matrix 
%CPSK_received = CPSK_trans.*fading + noise;  %nofading here!
CPSK_received = CPSK_trans + noise; 
%-------------------------- END OF TRANSMITTER --------------------------

%-------------------------- START OF RECEIVER ---------------------------


    CPSK_chaos_decoded = CPSK_received .* CPSK;

    accumulator = sum(CPSK_chaos_decoded);
    
    for ii = 1:2*N, 
        if accumulator(1,ii) < 0,
            output_vector(1,ii) = 0;  %unipolar
        else
            output_vector(1,ii) = 1;
        end;
    end;
        
  %%%Jimmy's decoding algorithm========================================
  tblen=15;      %the decoder depth is accounted by 5 times the constraint length(5*3=15)
  decoded=vitdec(output_vector,trel,tblen,'cont','hard');
  %====================================================================      
   code_output=convenc(decoded,trel);  %input vector binary 0 and 1     
        %compare3=xor(code_output(1,1+2*tblen:2*N),code(1,1:2*N-2*tblen));
        %compare4=xor(decoded(1,1+tblen:N),input_vector_uncoded(1,1:N-tblen));
%        
          jj = 1;
           for ii = 2*tblen+1:2*N,  
               
            %error calculation
            if code_output(1,ii) ~= code(1,jj),
                total_bit_errors(hh) = total_bit_errors(hh) + 1;
            else
                total_bit_errors(hh) = total_bit_errors(hh);
            end;
            jj = jj+1;
            end;        
      
%            jj=1;
%            for ii = tblen+1:N,  
%                
%             %error calculation
%             if decoded(1,ii) ~= input_vector_uncoded(1,jj),
%                 total_bit_errors(hh) = total_bit_errors(hh) + 1;
%             else
%                 total_bit_errors(hh) = total_bit_errors(hh);
%             end;
%             jj = jj+1;
%             end; 

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
    
    BER_eq(zz) = 0.5*erfc(( (psi/(SF/2) + ((users - 1)/(SF/2)) + ( SNR_lin(zz)./2).^-1)).^-.5);                                   %BER theory from CPSK 
    
end


semilogy(SNR_dB, BER_eq, 'rv-'); 
grid on;

legend('Simulated CPSK ','Theoretical CPSK - AWGN channel');

toc;

