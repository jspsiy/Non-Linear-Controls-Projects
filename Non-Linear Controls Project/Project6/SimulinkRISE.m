clear ;
clc
close all;

%% Parameters

%% GAINS

% k = 50;
% kn = 10;
% ks = 1;
% alpha = 10;
% gamma1 = 5;
% gamma2 = 10;
%k = 5;%[20,0;0,10];
%kn =1;
%ks = 1;
% beta1=0.2;
% alpha =0.1; 
% alpha2=alpha;
% gamma1 =500;
% gamma2 =50;
% wHatBar=5;
% vHatBar=5;
k =60;
kn=60;
ks=1;%not necessary ks=kn
beta1=0.2;
alpha=0.6;
alpha2=150;
gamma1=15000;
gamma2=30;
wHatBar=6;
vHatBar=6;
 v_initial =zeros(7,5);% 1-2*rand(7,5);
w_initial =zeros(6,2);%1-2*rand(6,2);
mew_initial = zeros(2,1);
f_hat_initial=zeros(2,1);
e0=[1;0];
eDot0=[0;0];
r0=eDot0+alpha*e0;



%% Run Simulation
[t,~,states,u,qd,f_hat,f,w_hat,v_hat] = sim('NNRISE');

%% Analysis/Plot

states = [states(:,3:4),states(:,1:2)]; % reorder states to [pos1, pos2, vel1, vel2]
error = qd-states;
eDot=error(:,3:4);
e=error(:,1:2);
r=eDot+alpha*e;
%rms = sqrt(mean(error.^2))*180/pi;
 
figure(1)
hold on;
plot(t,e,"-","LineWidth",2);
plot(t,r,":","LineWidth",2);
hold off;
title('tracking error')
legend('e1','e2','r1','r2')
xlabel('Time (s)')
ylabel('e and r')


figure(2)
hold on;
plot(t,u(:,1))
plot(t,u(:,2))
title('Link Controls')
legend('input1','input2')
xlabel('Time (s)')
ylabel('Control Torque (N-m)')
hold off;



figure(3)
title('F vs FHat')
subplot(1,2,1)
hold on;
plot(t,f(:,1),'-','LineWidth',2)
plot(t,f_hat(:,1),'--','LineWidth',2)
hold off;
legend('f1','f_h1');
subplot(1,2,2)
hold on;
plot(t,f(:,2),'LineWidth',2)
plot(t,f_hat(:,2),':','LineWidth',2)
hold off;
legend('f2','f_h2');

xlabel('Time (s)')
hold off;


figure(4)

[i,j,k]=size(w_hat);
hold on;
for x=1:i
    for y=1:j
        plot(t,squeeze(w_hat(x,y,:)),'--','LineWidth',2)
    end
end
title('What');
legend();
hold off



figure(5)
title('Vhat');
[i,j,k]=size(v_hat);
hold on;
for x=1:i
    for y=1:j
        plot(t,squeeze(v_hat(x,y,:)),'--','LineWidth',2)
    end
end
legend();
hold off




figure(6)
hold on;
plot(t,(f(:,1)-f_hat(:,1)))
plot(t,(f(:,2)-f_hat(:,2)))
hold off;
title('Ftilde')
legend('error1','error2')
xlabel('Time (s)')
ylabel('Percent Error (%)')








  

 

