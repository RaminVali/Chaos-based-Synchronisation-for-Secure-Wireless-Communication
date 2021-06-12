%%% 

tic;


clear;
clc;
%close all;

L = 1000; 
sequence_length = L   ; %chaotic sequence Length;




Rb = 1;  % bit rate = bits/sec
Rc = 1;   % chip rate = chips/bit   in other words spreading factor
Rs = 2;       % sample rate = samples/chip  

Rt = Rb*Rc*Rs;  %total rate



%What sort of sequence do you want correlated? Equal that to 1;
PN = 1;
CHAOS = 0;




%     %FOR CHAOS
%     x= zeros(1,L);                  % make the original x vector
%      initial_condition = 0.1133;
%      x(1) = initial_condition;
%      for ii = 2:length(x),
%         x(ii) = 1- 2*x(ii-1).^2;
%      end;



if PN ==1
    %FOR PN
    LFSR_connections = 5;
    sequence_type = 1;
    [x1] = m_sequence_generator(LFSR_connections, sequence_type);

    
    x = zeros(1,length(x1)*Rs);

   
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
    Rt = Rb*Rc*Rs;  %total rate
   time = 0:(1/Rt):(1/Rt)*(L*Rs-1);
   figure; stem(time,x,'r.');
end;




if CHAOS ==1
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
        %===========
end


%%% correlation
y_corr = [x,x,x];
R = zeros(2*length(x)-1,Rs);


index = 1;
while(index <= 3)  %index is for the delays, There are 3 indecies, 1 for lag, 2 for on time and 3 for half a chip advance. 
    y_ptr = index;
    ii = 1;
    while (y_ptr <(length(y_corr)-length(x))),
        R(ii,index) = sum (x(1:length(x)).*y_corr(y_ptr:y_ptr+length(x)-1));    
         y_ptr = y_ptr + 1; 
         ii = ii+1;
    end 
    index = index +1;
end
lags = (1:1:2*length(x)-1) - length(x); lags = lags.'; %%from Salim's code. This is the x axis which is scaled for the different sampling sizes. 
% the goal is to make the  final x axis look like the x axis of one sample
% perchip. but just simply have mopre samples in the space between them.

% lags =  1:length(R(:,1));
% lags = lags.';


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

%%%========================================================================
%%%===========================error calculation============================
%%%========================================================================
error = Rd - Ra;






%%%========================================================================
%%%=================================GRAPHICS===============================
%%%========================================================================


  % Create plot
%  plot(lags./(Rs),R(:,2)./max(R(:,2)),'DisplayName','The on time auto Correlaton function');
 % plot(lags./(Rs),R(:,2),'DisplayName','The on time auto Correlaton function');
%  hold on;

% Create plot
plot(lags./(Rs), Rd./max(Rd), 'DisplayName','Delayed','Color',[0 1 0]);
% plot(lags./(Rs), Rd, 'DisplayName','Delayed','Color',[0 1 0]);
hold on;
% Create plot
plot(lags./(Rs), Ra./max(Ra),'DisplayName','Advanced','Color',[1 0 0]);
% plot(lags./(Rs), Ra,'DisplayName','Advanced','Color',[1 0 0]);
hold on;
% plot(lags./(Rs), error./max(error),'DisplayName','Error','Color',[0 0 0]);
% plot(lags./(Rs), error,'DisplayName','Error','Color',[0 0 0]);
% Create legend
legend('show');
xlabel('\tau (in T_c)','fontsize',16);
ylabel('Normalized R(\tau)','fontsize',16);
%title('R(\tau), Erlay, Late and on time with 2 Samples/Chip- PN 127','fontsize',16);
axis([-10,10 -2 2]);

% plot(lags,R(:,2))
% hold on;
% plot(lags,error,'r')
% hold on;


%axis([-2 2 -1 1])
toc;

