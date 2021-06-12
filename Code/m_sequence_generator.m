function [sequence] = m_sequence_generator(no_of_taps, type)

if (no_of_taps == 2),
    taps = [1 1];
elseif (no_of_taps == 3),
    taps = [1 0 1];
elseif (no_of_taps == 4),
    taps = [1 0 0 1];
elseif (no_of_taps == 5),
    taps = [1 0 0 1 0; 1 1 1 1 0; 1 1 0 1 1];
elseif (no_of_taps == 6),
    taps = [1 0 0 0 0 1; 1 1 0 0 1 1; 1 1 0 1 1 0];
elseif (no_of_taps == 7),
    taps = [
        1 0 0 0 0 0 1; 1 0 0 0 1 0 0; 1 0 0 0 1 1 1;1 0 0 1 1 1 0; 
        1 1 0 1 0 1 0; 1 1 0 0 1 0 1; 1 1 1 0 0 1 0;1 1 1 1 0 1 1; 
        1 0 1 1 1 1 1];
elseif (no_of_taps == 8),
    taps = [
        1 0 0 0 1 1 1 0; 1 0 1 1 0 1 0 0; 1 0 1 1 0 0 1 0; 1 0 0 1 0 1 0 1; 
        1 0 1 1 0 0 0 1; 1 1 1 0 0 0 0 1; 1 1 1 1 0 0 1 1; 1 0 1 0 1 1 1 1];
elseif (no_of_taps == 9),
    taps = [
        1 0 0 0 0 1 0 0 0; 1 0 0 1 0 1 1 0 0; 1 1 0 0 1 1 0 0 0; 1 1 0 0 0 1 0 0 1;
        1 0 0 0 1 0 1 1 0; 1 1 0 1 1 0 0 0 0; 1 1 1 0 0 0 0 1 0; 1 0 0 1 1 1 0 1 1;
        1 0 1 1 0 1 1 0 1; 1 1 1 1 1 0 1 0 0];
elseif (no_of_taps == 10),
    taps = [
        1 0 0 0 0 0 0 1 0 0; 1 0 1 0 0 0 0 1 1 0; 1 0 0 0 0 0 1 1 0 1; 1 0 1 0 0 1 0 0 0 1;
        1 0 1 0 0 1 1 0 0 0; 1 1 0 0 0 0 1 0 0 1; 1 0 1 0 0 0 1 1 0 0; 1 0 0 0 0 1 0 1 1 0;
        1 0 0 0 0 1 0 0 1 1; 1 1 0 0 0 0 1 0 1 0];
elseif (no_of_taps == 11),
    taps = [
        1 0 0 0 0 0 0 0 0 0 1; 1 0 0 1 0 0 1 0 0 1 0; 1 0 0 0 1 0 0 0 1 1 0;
        1 0 0 0 0 0 1 0 1 1 0; 1 1 0 0 0 0 0 0 1 1 0; 1 0 0 0 0 1 1 0 0 0 1;
        1 0 0 0 0 0 1 0 1 0 1; 1 0 1 0 0 0 0 1 0 0 1; 1 0 0 1 0 1 0 0 0 1 0;
        1 0 1 1 0 0 0 0 1 0 0];
elseif (no_of_taps == 12),
    taps = [
        1 0 0 0 0 0 1 0 1 0 0 1; 1 0 0 1 0 0 0 0 0 1 1 0; 1 1 1 0 0 0 0 1 0 0 1 1;
        1 1 0 0 0 0 1 0 1 0 1 1; 1 1 0 1 0 1 1 1 0 0 0 0; 1 1 0 1 0 0 0 1 0 1 0 1;
        1 1 0 1 1 1 0 0 1 0 0 0; 1 1 0 1 0 1 1 1 0 0 0 0; 1 0 0 1 1 0 0 0 0 1 1 1;
        1 0 1 1 1 0 1 0 0 0 1 0];
elseif (no_of_taps == 13),
    taps = [
        1 0 0 0 0 0 0 0 0 1 1 0 1; 1 0 0 1 1 0 1 0 1 1 0 0 0; 1 0 1 0 0 1 1 0 0 1 0 0 1;
        1 1 0 0 0 1 1 1 1 0 0 0 0; 1 0 0 0 1 1 1 0 1 0 0 0 1; 1 1 0 0 0 0 0 1 1 1 1 0 0;
        1 1 1 0 1 0 0 0 1 0 1 0 0; 1 1 1 0 0 0 0 0 1 0 0 1 1; 1 1 0 0 1 1 0 0 0 1 0 1 0;
        1 0 0 0 0 1 1 0 0 1 1 1 0];
else
    sequence = 0;
    return;
end;

m = size(taps,1);   %Number of rows.
rand('state',sum(100*clock));
if (m > 1),
    x = round((m - 1)*rand(1,1) + 1);   %Randomly select taps combination.
    taps = taps(x,:);
end;

if (type == 0),
    sequence = ss_mlsrs(taps);      %Unipolar sequence ... {0,1}
else
    sequence = 2.*ss_mlsrs(taps)-1; %Bipolar sequence ... {-1,+1}
end;


%=====================================================
%           GENERATION OF M-SEQUENCES 
%=====================================================

%Function obtained from (page 408):
%Proakis, Salehi, and Bauch, Contemporary Communication Systems Using MATLAB and Simulink, 
%2nd ed. Toronto, Ontario: Thomson - Brooks/Cole, 2003.

function [seq] = ss_mlsrs(connections)

m = length(connections);
L = 2^m - 1;
registers = [zeros(1,m-1) 1];
seq(1) = registers(m);

for i=2:L,
    new_reg_cont(1) = connections(1)*seq(i-1);
    for j=2:m,
        new_reg_cont(j)= xor(registers(j-1),connections(j)*seq(i-1));
    end;
    registers = new_reg_cont;
    seq(i) = registers(m);
end;