%function twoLinkRobotAdaptive
close all
clear all
%Set up parameters for sim
p1       = 3.473;
p2       = 0.196;
p3       = 0.242;
f1       = 5.3;
f2       = 1.1;
global tau tilde eps e20 ar;
tau=[];
tilde=[];
eps=[];
% Stacked parameter vector
theta    = [p1;p2;p3;f1;f2];

% Simulation final time
tf   = 60;

% Initial condition vector (X0 must be same size and "form" as X and XDot below)
% (i.e., in this sim, X0 = [e0;r0;thetahat0])
% X0   = [4;10;%e1 @2
%         0;0;%e2  @4
%         3;2; %r @6
%         1;1;1;1;1;%thetahat @11
%         0;0;0;0;0; %omegarow1 @16
%         0;0;0;0;0;%omegarow2 @21
%         0;0;0;0;0;%Wrow1@26
%         0;0;0;0;0;%Wrow2 @31
%         0;0;  % ufd states @33
%         0;0]; % mew states @35
    
    X0   = [4;10;%e1 @2
        3;2;%e2  @4
        0;0; %mew @6
        1;1;1;1;1;%thetahat @11
        ]; % mew states @35
e20=X0(3:4);
% Options for integration function
opts = odeset('RelTol',1e-3,'AbsTol',1e-3);

qd0 = [cos(0.5*0);2*cos(0)];
qdDot0 = [-0.5*sin(0.5*0);-2*sin(0)];
e0=[X0(1);X0(2)];
r0=[X0(3);X0(4)];
q0=e0+qd0;
a=1.5;
qDot0=r0+qdDot0-a*e0;
c20=cos(q0(2));
c2d0=cos(qd0(2));

M0qDotY= [q0(1), q0(2)     ,2*c20*q0(1)+c20*q0(2),0,0;
            0  ,q0(1)+q0(2),c20*q0(1)            ,0,0];
Md0qdDotY= [qd0(1), qd0(2)       ,2*c2d0*qd0(1)+c2d0*qd0(2),0,0;
            0     , qd0(1)+qd0(2),c2d0*qd0(1)              ,0,0];

% Integrate (you can send the paramters theta to the dynamics as seen below)
[t,STATES] = ode45(@(t,X) twoLinkdynamics(t,X,M0qDotY,Md0qdDotY,theta),[0 tf],X0,opts);

% Set up desired trajectory data for plots (enter desired trajectory for your simulation)
qd = [cos(0.5*t) 2*cos(t)]';

% Parse integrated states (STATES is the same "form" as X0)
% (i.e., in this sim, STATES = [e r thetahat] over all time);
e1  = STATES(:,1:2)';
e2  = STATES(:,3:4)';
thetaHat = STATES(:,7:11)';
%r  = STATES(:,5:6)';
%thetaHat = STATES(:,7:11)';

% Compute x from e and xd for plotting purposes
q  = qd-e1;

% Plot the actual vs desired trajectories

% Plot the actual vs desired trajectories
subplot(2,2,1)
plot(t,qd,'-','LineWidth',2)
hold on
ax = gca;
ax.ColorOrderIndex = 1;
plot(t,q,':','LineWidth',2)
title('q vs qd');
hold off

% Plot the filtered tracking error
subplot(2,2,2)
hold on

plot(t,e1,':','LineWidth',2)
plot(t,e2,':','LineWidth',2)
legend("e1a","e1b","e2a","e2b");
title('Error')
hold off

% Plot the adaptive estimates vs actual parameters
 subplot(2,2,3)
 plot(t,repmat(theta,1,length(t)),'-','LineWidth',2)
 hold on
 ax = gca;
 ax.ColorOrderIndex = 1;
 plot(t,thetaHat,':','LineWidth',2)
 title('Thetahats');
 hold off
% 
% subplot(2,2,3)
% plot(ar(:,1),ar(:,2),'--','LineWidth',2)
% hold on
% plot(ar(:,1),ar(:,3),'-','LineWidth',2)
% legend("r-error1","r-error2");
% title('r-error')
% hold off


subplot(2,2,4)
plot(tau(:,1),tau(:,2),'--','LineWidth',2)

hold on
plot(tau(:,1),tau(:,3),'-','LineWidth',2)
legend("torque1","torque2");
title('input')
hold off

 figure(2)
 hold on
 plot(tilde(:,1),tilde(:,2))
 plot(tilde(:,1),tilde(:,3))
 plot(tilde(:,1),tilde(:,4))
 plot(tilde(:,1),tilde(:,5))
 plot(tilde(:,1),tilde(:,6))
 legend('tilde1','tilde2','tilde3','tilde4','tilde5')
 title("thetatilde")
% figure(3)
% hold on
% plot(eps(:,1),eps(:,2))
% plot(eps(:,1),eps(:,3))
% legend('E1','E2')
% title("Epsilon")

function [XDot] = twoLinkdynamics(t,X,M0qDotY,Md0qdDotY,theta)
global tau tilde eps e20 ar;
% Parse parameter vector
p1 = theta(1);
p2 = theta(2);
p3 = theta(3);
f1 = theta(4);
f2 = theta(5);

% Select gains for controller
K        = 50*eye(2);%[15,0;0,1.5];%9; %Enter a number
a1        = 2;%0.15; %Enter a number
a2        =2;

beta=0.5;

% Desired trajectory and needed derivatives
qd       = [cos(0.5*t);2*cos(t)];
qdDot    = [-0.5*sin(0.5*t);-2*sin(t)]; %Enter the expression
qdDotDot = [-0.25*cos(0.5*t);-2*cos(t)]; %Enter the expression


% Parse current states (X is the same size and "form" as X0)
% (i.e., in this sim, X = [e;r;thetahat])
e1       = [X(1);X(2)];
e2       = [X(3);X(4)];
mew=[X(5);X(6)];
%thetaHat = [X(7);X(8);X(9);X(10);X(11)];
%ohmyd=[X(12:16)';X(17:21)'];%omega parameter
%wyd=[X(22:26)';X(27:31)'];
%uf=[X(32);X(33)];

% Compute current x and xDot for convenience
q        = qd-e1;
e1Dot    =e2-a1*e1;
%e2Dot    = r - a2*e2 ;
qDot=qdDot+a1*e1-e2;
%qDotDot=qdDotDot-r+a2*e2+a1*e1Dot;
%e1DotDot=r-a1*e1Dot-a2*e2;
% Compute cos(x2) and sin(x2) for convenience
c2       =cos(q(2));
s2       =sin(q(2));
c2d      =cos(qd(2));
s2d      =sin(qd(2));

% Compute current matrices for the dynamics
M        = [p1+2*p3*c2 p2+p3*c2;p2+p3*c2 p2];
%mDot=[-2*p3*s2*qDot(2),-p3*s2*qDot(2);-p3*s2*qDot(2),0];
Vm       = [-p3*s2*qDot(2) -p3*s2*(qDot(1)+qDot(2));p3*s2*qDot(1) 0];
%VmDot=[-p3*c2*qDot(2)^2-p3*s2*qDotDot(2),-p3*c2*qDot(2)*(qDot(1)+qDot(2))-p3*s2*(qDotDot(1)+qDotDot(2));p3*c2*qDot(2)*qDot(1)+p3*s2*qDotDot(1),0];
fd       = [f1 0;0 f2];
tauD=[0.5*cos(0.5*t);sin(t)];
%tauDDot=[-0.25*sin(0.5*t);cos(t)];



%u=mew;
%e2Dot=M\(-Vm*qDot-fd*qDot+u-tauD)-qdDotDot+a1*e1;
%muDot=(K+1)*(e2Dot+a2*e2)+beta*sign(e2);
%(K+1)*a1*e2 - beta*sign(e2);
u        =(K+1)*e2-(K+1)*e20+mew;%((K+1)*e2 - (K+1)*e20 + mew);%+Y*thetaHat;%-K*r; % -e Enter the expr..ession
e2Dot=qdDotDot-M\(-Vm*qDot-fd*qDot-tauD+u)+a1*e1Dot+a2*e2;
r=e2Dot+a2*e2;
muDot=(K+1)*a2*e2 + beta*sign(e2);
thetahatDot=zeros(5,1);

%S=M*qdDotDot+Vm*qDot+fd*qDot+M*(a1*e1Dot+a2*e2)-Yd*theta;
%SDot=mDot*qdDotDot+M*qdDotDotDot+VmDot*qDot+Vm*qDotDot+fd*qDotDot+mDot*(a1*e1Dot+a2*e2)+M*(a1*e1DotDot+a2*e2Dot)-YdDot*theta+Yd*thetaHatDot;
%r=M\(Yd*theta+S+tauD-u);
%rDot=inv(M)*(-mDot*r+YdDot*theta-Yd*thetaHatDot+SDot+tauDDot-tauDot)
thetaHat       = [X(7);X(8);X(9);X(10);X(11)];
tilda=theta-thetaHat;
tau=[tau;t(length(t)),u'];
ar=[ar;t(length(t)),r'];
tilde=[tilde;t(length(t)),tilda'];
%eps=[eps;t(length(t)),Epsilon'];
% Compute current closed-loop dynamics

%rDot        = M\(-Vm*qDot-fd*qDot-TauD+u)-qdDotDot+a1*e1Dot+a2*e2; %Enter the expression
%rDot        = qdDotDot-M\(-Vm*qDot-fd*qDot-TauD+u)+a1*e1Dot+a2*e2;
%Ntilde=-Yd*gamma*transpose(YdDot)*r;
%rDot=M\(-1/2*mDot*r+YdDot*thetatilda+Ntilde+Nd-muDot-e2);
% eye(5)*transpose(Yd)*r;%Enter the expression


% Stacked dynamics vector (XDot is the same size and "form" as X)
%XDot        = [e1Dot;e2Dot;rDot;thetaHatDot;ohmydMatDot;wydMatDot;ufdDot;muDot];
XDot        = [e1Dot;e2Dot;muDot;thetahatDot;];
end