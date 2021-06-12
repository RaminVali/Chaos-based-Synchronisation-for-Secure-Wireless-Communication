%%%DLL_v1.m   Date:10/12/2007  By RV
%Description:

clear;
clc;
%Parameters

Rb = 200;  % bit rate = bits/sec
Rc = 100;   % chip rate = chips/bit   in other words spreading factor
Rs = 1;      % sample rate = samples/chip  

Rt = Rb*Rc*Rs;  %total rate
pilot_period= 1000*Rs;  %We want the pilot and the correlation lengths in samples
corr_period = 200*Rs;


offset = 100;           %in samples






% make the original x vector
x = zeros(1,pilot_period);                
 initial_condition = 0.1133;
 x(1) = initial_condition;
 for ii = 2:length(x),
    x(ii) = 1 - 2*x(ii-1).^2;
    
 end;
 
 for ii = 2:length(x),
   
    if   x(ii)>= 0
      
        x(ii) = 1;
    else
        x(ii) = -1;
    end;
 end;
 x(1) = 1;
 
 
 
 %interpolate the x vector
 %x_interpolated = interp(x,Rs);  %% generate interpolated signal
 
 %assign tx-array and rx_array
 tx_array = x;
 rx_array = tx_array;
 

 tx_ptr = offset;
 tx_ptr_s = offset; 
tx_ptr_f = offset + corr_period; 
 
 rx_ptr = 0;
 rx_ptr_s =  200; %starting rx pointer
 rx_ptr_f = rx_ptr_s + corr_period; 
 %finishing rx pointer
  if rx_ptr_f >= pilot_period
             rx_ptr_f = mod(rx_ptr_f,pilot_period) ;
  end;







%------------------ Find the Auto-correlation of the pilot ----------------

for jj = 1:pilot_period,
    % Take the next T samples:
    for ii = 1:corr_period,
        % Modulo the long way since matlab has no concept of A(0)
        if rx_ptr >= pilot_period,
            rx_ptr = 0;
        end;
        if tx_ptr >= pilot_period,
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

    if rx_ptr >= pilot_period,
         rx_ptr = 0;
    else 
         rx_ptr = rx_ptr + 1;
    end;
end;



plot(R,'r');


 
 
