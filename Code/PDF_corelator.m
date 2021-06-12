
clc;

clf;


ii = 0;
users = 6;
%result = zeros(1,3000*users);
load Probability.mat
x = prob(1,500:2501);
 
 bins = linspace (-1,1,length(x));
 plot(bins,x);
 title(['Logistic map PDF of ', int2str(1), ' users ' ])
 xlabel('Chaotic Value');
 ylabel('Probability');

  result = conv(x,x); 
    bins = linspace(-2,2,length(result));
  
  figure;plot(bins,result);
title(['Logistic map PDF of ', int2str(2), ' users'])
xlabel('Chaotic Value');
 ylabel('Probability');

 
for ii = 3:users,
    
  result = conv(x,result);
  bins = linspace(-ii,ii,length(result));
 figure; plot(bins,result); title([' Logistic map PDF of ', int2str(ii), ' users']);
 xlabel('Chaotic Value');
 ylabel('Probability');
end