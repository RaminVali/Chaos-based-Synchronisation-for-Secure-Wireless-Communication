clear;
clc;

Length = 1000000;

SNR_dB = 0;
tx_array = zeros(1,Length);

sampling_freq = 1000;
%SF = [10,100,1000];  % sending 10 bits now (SF = 100)

sf = 100;


chaos = zeros(1,Length);

chaos(1,1) = 0.231;
for ii = 2 :length(chaos),
chaos(1,ii) = 1-2.*(chaos(ii-1).^2);
end;



 % Calculate the noise - do the long way
   
Pavg = 0.5;%sum(tx_array.^2) / (length(tx_array));   % Variance of the logistic map = 0.5
Pavg = sum(chaos.^2) / (length(tx_array)); 
    % Eb is the total energy of the signal divided by the number of bits:
    Etot = Pavg .* (length(tx_array)/sampling_freq);
    Eb = Etot./(length(tx_array)./sf);
    % and Ec is the total energy of the signal divided by the number of chips
    Ec = Etot./(length(tx_array));
    SNR_lin = 10.^(SNR_dB/10);
    No = Eb/SNR_lin;
    sigma_squared = No*sampling_freq/2;  %%Fs/2 is the BW
    

Etot_2 = sum(chaos.^2)



% Eb = zeros(1,1000);
% Eb(1,1) = sum(chaos(1,SF).^2);
% 
% for(kk =1:length(SF)),
%     
%     for jj = 2:length(Eb),
%         Eb(1,jj) = sum(chaos(1,((jj-1)*SF(kk)):jj*SF(kk)).^2);
%     end;
%     plot(Eb);
%     plot(Eb./SF(kk),'r');
%     hold on;
% end;
