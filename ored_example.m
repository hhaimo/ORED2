tstep=1e-3; % time discretization step
tf=3;
t=[0:tstep:tf]; % time instants

L=1e3; % Bound on 3rd-order derivative of f
N=1e-1; % Noise bound (unknown to the differentiator)

R2=100; %130; % Initial second derivative
R1=0; %100; % Initial first derivative
f=L*(t.^3)/6 + R2/2*(t.^2) + R1*t;
derf1=L*(t.^2)/2 + R2*t + R1; 
derf2=L*t + R2; 

% Lowest worst-case bounds
M1=3*18^(1/3)/2*L^(1/3)*N^(2/3);
M2=2*12^(1/3)*L^(2/3)*N^(1/3);

% Convergence time
tconv=(180*N/L)^(1/3);

% noise and input signal
seed=815; % Random Number Generator seed
rng(seed);
% Noise realization
eta=N*ones(size(t)).*(t<1)+ ...
    (N-(L*2.1)/6*(t-1).^3).*(t>=1).*(t<1.2)+ ...
    (rand(size(t))-.5)*2*N.*(t>=1.2).*(t<2) + ...
    N*sin(500*t).*(t>=2);
eta=max([-N*ones(size(t));eta]); % saturate below
eta=min([N*ones(size(t));eta]);  % saturate above
u=f+eta; % Measurements

%% ORED 
%
Nest=0;
clear Thats Ths N1hats N2hats Nhats y1 y2
Nhat_prev=[0 0];
for ii=4:1:length(u) % Index t
    Nest=Nest+1;
    [N1hat,N2hat]=noise_estimate(u(1:ii),tstep,L);
    N1hats(Nest)=N1hat;
    N2hats(Nest)=N2hat;
    N12hat=max([N1hat N2hat]);
    Nhat=max([N12hat Nhat_prev]);
    Nhats(Nest)=Nhat;
    Nhat_prev=[Nhat_prev(2) N12hat];
    That1=(180*Nhat/L).^(1/3);
    That=min([(ii-1)*tstep,That1]);
    Thats(Nest)=That;
    % That may be not a multiple of 3*tstep, correction:
    Th3_ind=max([floor(That/(3*tstep)) 1]); % integer index of That/3, not less than 1
    Th_ind=Th3_ind*3; % integer index of That
    Th=Th_ind*tstep;  % corrected That
    Ths(Nest)=Th;
    y1(Nest)=(8*u(ii)-9*u(ii-Th3_ind)+u(ii-Th_ind))/(2*Th);
    y2(Nest)=(6*u(ii)-9*u(ii-Th3_ind)+3*u(ii-Th_ind))/Th^2;
end

ored_time=t(4:1:end)';
ored_de1=abs(derf1(4:1:end)-y1)'; % ORED absolute error on first derivative
ored_de2=abs(derf2(4:1:end)-y2)'; % ORED absolute error on second derivative

%% Plots

figure(1)
clf
plot(t,u), hold on
plot(t,f)
title('Input and signal')
grid on

figure(2)
clf
plot(t(4:1:end),abs(eta(4:1:end))), hold on
plot(t(4:1:end),Nhats)
title('Noise absolute value and estimate')
grid on

figure(3)
clf
plot(t(4:1:end),Thats),hold on
plot(t(4:1:end),Ths)
title('T hat and corrected T hat')
grid on

figure(4)
clf
plot(t(4:1:end),y1,'b'), hold on
plot(t,derf1,'k')
plot(t,derf1+M1,'k--')
plot(t,derf1-M1,'k--')
%plot(y_ired.Time,y_ired.Data(:,1),'r')
cur_axis=gca;
plot([tconv tconv],cur_axis.YLim,'k--')
title('first derivative')
grid on

figure(5)
clf
plot(t(4:1:end),y2,'b'), hold on
plot(t,derf2,'k')
plot(t,derf2+M2,'k--')
plot(t,derf2-M2,'k--')
%plot(y_ired.Time,y_ired.Data(:,2),'r')
cur_axis=gca;
plot([tconv tconv],cur_axis.YLim,'k--')
title('second derivative')
grid on

figure(6)
clf
plot(t(4:1:end),ored_de1,'b') %abs(derf1(4:1:end)-y1),'b')
hold on
%plot(t,ired_de1,'r') %abs(derf1'-y_ired.Data(:,1)),'r')
plot([0 tf],[M1 M1],'k--')
cur_axis=gca;
plot([tconv tconv],cur_axis.YLim,'k--')
title('First derivative absolute error')
grid on

figure(7)
clf
plot(t(4:1:end),ored_de2,'b') %abs(derf2(4:1:end)-y2),'b')
hold on
%plot(t,ired_de2,'r') %abs(derf2'-y_ired.Data(:,2)),'r')
plot([0 tf],[M2 M2],'k--')
cur_axis=gca;
plot([tconv tconv],cur_axis.YLim,'k--')
title('Second derivative absolute error')
grid on

figure(8)
clf
plot(t,eta)
title('Noise')
grid on

