%function ZOHSample, samples the \input signal ZOH fashion with the given
%number of samples per chip.

function [x_ZOH] = ZOHSample(x,Rs,sequence)
clc;
if sequence == 1
initial_condition = 0.1133;

 %% for the zero order hold sampling (square sampling method)
        x_zoh = zeros(Rs,length(x)-1); % to get and exact sequence length
                                              % we want, matlab arrays start at 1
        x_zoh(:,1)=initial_condition;   %initial condition
        
        for (ii = 1:sequence_length-1),
            x_zoh(:,ii+1) = 1- 2.* (x_zoh(:,ii).^2); %makes a matrix of dimentions(sample,sequencelength)  the coloums contain the same numbers
        end;                                 
         
        y_zoh  = [x_zoh(:,1).'];        % get the first colum and make it into the first row
           
        
        for (jj  = 2:sequence_length),      
           y_zoh = [y_zoh,x_zoh(:,jj).'];      %keep appending to this row to make the transmission vector.      
        end;
        x_ZOH = y_zoh;
        
        
elseif sequence == 0,
    x_ZOH = zeros(1,length(x)*Rs); % Introduce ZOH sampling of the signal 
        jj = 1;
        for mm = 1:length(x)
           for nn = 1:Rs
               x_ZOH(jj) = x(mm);
               jj = jj+1;
           end
            nn = 1; 
        end;
    
end;