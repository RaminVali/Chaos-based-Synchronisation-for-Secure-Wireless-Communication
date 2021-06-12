
%21/6/07
%%Ramin Vali
%%%This program takes two arays ( chaotic one already generated) and gives
%%%the correlation ( cross or auto) results and plots them. The workings
%%%are compared to the matlab xcorr function and the error between them is
%%%0 for all the purposes needed.

%%for the matters of calculaing  the theoretical auto correlation, we need
%%the secondary  array to be padded with the first array for 3 times,
%%becayse these functions are perodic rather than causal( say a PN sequence
%%with N = 7. so we buffer up with 3 of thses and the auto correlation
%%function will come up like haykin and other erefernces. We had the same
%%problem with chaos. When we slide them past each other, we need to have
%%more of one. hence the point of the circular buffer,

%%%try making a circular buffer.




clc;
clear;
load 'store';

%%%===================~~~~~~~~~~~~~~value generation ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 x = zeros(1,100);
% x1 = zeros(1,1000);
R = zeros(1,2*length(x)-1);
 

ii = 1;
y_ptr = 2;

 x(1) = 0.2331;
 for kk = 2:length(x),
     x(kk) = 1- 2*x(kk-1).^2;
%       if x(kk)<=0      %%this is used when unipolar is desired 
%           x(kk) = -x(kk);
%      end;
 end;
%x = round(store);

% LFSR_connections = 13;
% sequence_type = 0;
% [x] = m_sequence_generator(LFSR_connections, sequence_type);
%  

%   x1(1) = 0.231;
%  for kk = 2:length(x),
%      x1(kk) = 1- 2*x1(kk-1).^2;
% %      if x(kk)<=0.5      %%this is used when unipolar is desired 
% %          x(kk) = -x(kk);
% %      end;
%  end;
%  
 %%%padding y witch zeros so we can make the algorithm easier The padding
 %%%size is the size of the secondary array both on the left and right of
 %%%it.(or more just to see the periodicity) ( look at hakin theory)
 %%%chapter 7
 y = [x,x,x];  %% padding with the array itself
%y = [zeros(1,length(x)) x zeros(1,length(x))];  %%padding with zeros
 %% buffer up with the periodic!!( signal again)

%%the main algorithm, basically while the pointer on the y array is less
%%than length of the larger array - length of the smaller array, the
%%correlation is calculated by multiplying all the members of the x array
%%member by member by the y array from the pinter to pointer +length of x.
%%This is done member by member.  then increment the pointer and start
%%again 
  while (y_ptr <(length(y)-length(x))),
    R(ii) = sum (x(1:length(x)).*y(y_ptr:y_ptr+length(x)-1));    
     y_ptr = y_ptr + 1; 
     ii = ii+1;
  end  
        
    
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~graphics~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


% Create figure
figure1 = figure;

% Create subplot
subplot1 = subplot(3,1,1,'Parent',figure1);
stem(x,'Parent',subplot1);
xlabel('Index','fontsize',16);
ylabel('Value','fontsize',16);
title('Chaotic Bipolar','fontsize',16);

% Create subplot
subplot2 = subplot(3,1,2,'Parent',figure1);
plot(R./max(R),'Parent',subplot2);
xlabel('Index','fontsize',16);
ylabel('Value','fontsize',16);
title('Autocorrelation Function','fontsize',16);
legend('My  algorithm');



% Create subplot
subplot3 = subplot(3,1,3,'Parent',figure1);
plot(xcorr(x),'r','Parent',subplot3);
xlabel('Index','fontsize',16);
ylabel('Value','fontsize',16);
title('Autocorrelation Function','fontsize',16);
legend('Matlab xcorr algorithm');







