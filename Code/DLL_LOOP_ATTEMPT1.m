%%% DLL LOOP ATTEMPT1       15/2/08

%modified 18/2/08  Increased the sample resultion Discovered the
%                   relationship between the accumulative error and the
%                   resolution of the trackingf loop S curve.
%modified 20/2/08   Discovered the relationship between the static error
%                   and the initial condition of the time estimate
clc;
clear;
%load TimingError.mat;
% Rs = 10;
% TimingError = TimingError.'./Rs;

 t = linspace(-0.5,0.5,1e4);


    TimingError = 1024.*t + 200*-0.5*randn(1,length(t));



N = 511;
Ap = 1;
m = 1024; % the slope ( for this case ONLY)
Td = zeros(1,40);
Td_hat = (zeros(1,length(Td)-1));

Td_hat(1) = 0;

x = 0.5 - rand;
      Td(1) = x;


for kk = 1:length(Td)
    
    if kk >1
        x = 0.5 - rand;
       
        Td(kk) = x + Td(kk-1);
        Td_pure(kk) = x;
     
    end 
    
   for zz = 1:length(t),
       if x >= t(zz) && x<t(zz+1),
           
          y = zz;
           break;
       end;
       
   end;

    Td_hat(kk+1) = Td_hat(kk)-((TimingError(y))./(Ap.^2.*(2*N+2)));
    
    
temp(kk+1) =  (((TimingError(y))./(Ap.^2.*(2*N+2))));
end;
% 
% Create figure
figure1 = figure('PaperSize',[20.98 29.68]);

% Create subplot
subplot1 = subplot(2,1,1,'Parent',figure1,'FontSize',18,'FontName','Times New Roman');
box('on');
hold('all');

Td_hat =  (Td_hat);
 plot(-Td_hat,'r--');
 hold on;
 plot(Td);
%  
xlabel('$$t/T_c$$','interpreter', 'latex','fontsize',30');
ylabel('$$\hat{T}_{d} \hspace{0.01cm}\textrm{,}\hspace{0.2cm}T_d$$','interpreter', 'latex','fontsize',30);
%title('$$\textrm{Cumulative}\hspace{0.2cm}\hat{T}_{d}  \hspace{0.2cm}\textrm{and}\hspace{0.2cm}T_d\hspace{0.2cm}\textrm{for each chip period}$$','interpreter', 'latex','fontsize',18);
 h = legend('$$\hat{T}_{d}$$','$$T_{d}$$');
 set(h, 'interpreter', 'latex','fontsize', 22);
 

% axis square

subplot2 = subplot(2,1,2,'Parent',figure1,'FontSize',18,...
    'FontName','Times New Roman');
box('on');
hold('all');

plot(temp,'r--');
hold on;
plot(Td_pure);

%axis square
%title('Non-cumulative T_d and T^_d for each chip period');
xlabel('$$t/T_c$$','interpreter', 'latex','fontsize',25');
ylabel('$$\hat{T}_{d} \hspace{0.01cm}\textrm{,}\hspace{0.2cm}T_d$$','interpreter', 'latex','fontsize',30);
%title('$$\textrm{Non-cumulative}\hspace{0.2cm}\hat{T}_{d}  \hspace{0.2cm}\textrm{and}\hspace{0.2cm}T_d\hspace{0.2cm}\textrm{for each chip period}$$','interpreter', 'latex','fontsize',18);
 h = legend('$$\hat{T}_{d}$$','$$T_{d}$$');
 set(h, 'interpreter', 'latex','fontsize', 22);
 
%xlabel('$$\frac{t}{T_c}$$','interpreter', 'latex','fontsize',25');

%  figure;
%    error = Td - Td_hat(1,1:length(Td));
%  plot(error,'k');
% grid minor;
