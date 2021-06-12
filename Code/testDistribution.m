% M-file to test the uniform distribution, to be used in diagnostics procedures.


clear;
%clc;
jj = 1;
kk = 1;
fid = fopen('Debug\uniform_distribution_BitShift.txt','r');
fgets(fid);
t = zeros(200, 2);
for i = 1:200
     t(i,:)= fscanf(fid, '%g',1);
end;
%t = fscanf(fid,'%g%g');
fclose(fid)

for ii = 1:200,
     if rem(ii,2) == 0
        a(jj) = t(ii);
        jj= jj+1;
     else
        b(kk) = t(ii);
        kk = kk+1;
     end;
end;

a = a.';
%b = b.';
b = 0.5:1:100;  % distribution measured at the centre
prob = a/sum(a);

bar(b/10, prob/max(prob), 1);
%plot(b/10, prob, 'r');
hold on;

