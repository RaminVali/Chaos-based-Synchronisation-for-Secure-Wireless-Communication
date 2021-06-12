%function that generates an array of logisitic map chaotic values

%Input  : the length of the vector needed.
%Output : the chaotic vector


function output = Chaos_gen(length)

IC =  0.123; %initial condition for the logictic map

output = zeros (length);
output(1) = IC;

for ii = 2:length,
    
    output(ii) = 1- 2*output(ii-1).^2;
 end;
 
 output = cast(output, 'double');