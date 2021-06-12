%clear;
clc;
figure;
MAX_DATAPOINTS = 10000;
% MAX_DATAPOINTS = input('please type the number of dataponits  ');
 

% for Rayleigh the negatives do not contribute.


kk = 1;
mm = 1;
pos_bin = 0;
neg_bin = 0;




% File to read in and plot the Rayleigh distribution

%fid = fopen('Debug\1Billion Sample runs\AWGN_BitShift_Results.txt', 'r');
 %fid = fopen('Fading_nofilter_10000.txt','r');
 
 %fid = fopen('fading_results_flitered1M.txt','r');  %% activate the
 %smaller bins to get the Rayleigh distribution for the filtered results
 
 %fid = fopen('fading_results.txt','r');
 %fid = fopen('gaussian_I.txt','r');
    %fid = fopen('gaussian_I.txt','r');


    
    chaos = zeros(1,100);

chaos(1,1) = 0.231;
for ii = 2 :length(chaos),
chaos(1,ii) = 1-2.*(chaos(ii-1).^2);
end;
    
    
    
    
%Data = fscanf(fid, '%f');
%Data = chaos;


%load('error');
% finding mean and variance of the data

%Data = randn(1,MAX_DATAPOINTS);

variance = var(error);
std_dev = sqrt (variance);
average = mean(error);


error = error.';
error = error./max(error);
%bins = 0:0.01:0.5;  % for Rayleigh ( filtered )

bins = -1.5:.0001:1.5;  %For Rayleigh non filtered or anything gAussian Just mind where the bins start 

% finds the number of the values betwean one value and the next, this is
% then put into the prob array to get the bins right.

for ii = 1:(length(bins)-1),
    prob(ii) = length(find(error>bins(ii)))-length(find(error>bins(ii+1)));
end;

bar(bins (1:(length(bins)-1)),prob/max(prob),'r');
%plot(Data)

hold on;







