%function twoLinkRobotAdaptive
close all
clear all
%Set up parameters for sim
p1       = 3.473;
p2       = 0.196;
p3       = 0.242;
f1       = 5.3;
f2       = 1.1;
global tau tilde eps;
tau=[];
tilde=[];
eps=[];
% Stacked parameter vector
theta    = [p1;p2;p3;f1;f2];

% Simulation final time
tf   = 60;

% Initial condition vector (X0 must be same size and "form" as X and XDot below)
% (i.e., in this sim, X0 = [e0;r0;thetahat0])
X0   = [4;10; %e
        3;2;  %r
        1;1;1;1;1;%thetahat @9
        0;0;0;0;0; %omegarow1 @14
        0;0;0;0;0;%omegarow2 @19
        0;0;0;0;0;%Wrow1@24
        0;0;0;0;0;%Wrow2 @29
        0;0;0;0;0; %desired omegarow1 @34
        0;0;0;0;0;% desired omegarow2 @39
        0;0;0;0;0;% desired Wrow1@44
        0;0;0;0;0;%desired Wrow2 @49
        0;0]; % uf states

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
e  = STATES(:,1:2)';
r  = STATES(:,3:4)';
thetaHat = STATES(:,5:9)';

% Compute x from e and xd for plotting purposes
q  = e + qd;

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
plot(t,r,'--','LineWidth',2)
plot(t,e,':','LineWidth',2)
legend("r1","r2","e1","e2");
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
figure(3)
hold on
plot(eps(:,1),eps(:,2))
plot(eps(:,1),eps(:,3))
legend('E1','E2')
title("Epsilon")

function [XDot] = twoLinkdynamics(t,X,M0qDotY,Md0qdDotY,theta)
global tau tilde eps;
% Parse parameter vector
p1 = theta(1);
p2 = theta(2);
p3 = theta(3);
f1 = theta(4);
f2 = theta(5);
ohmy=[X(10:14)';X(15:19)'];%omega parameter
wy=[X(20:24)';X(25:29)'];
ohmyd=[X(30:34)';X(35:39)'];
wyd=[X(40:44)';X(45:49)'];
% Select gains for controller
K        = [15,0;0,1.5];%9; %Enter a number
a        = .2;%0.15; %Enter a number
gamma    = zeros(5);%75;
gamma(1,1)    = 20;%75;
gamma(2,2)    = 10;%75;
gamma(3,3)    = 10;%75;
gamma(4,4)    = 60;%75;
gamma(5,5)    = 10;%75;
beta=0.5;
uf=[X(50);X(51)];
% Desired trajectory and needed derivatives
qd       = [cos(0.5*t);2*cos(t)];
qdDot    = [-0.5*sin(0.5*t);-2*sin(t)]; %Enter the expression
qdDotDot = [-0.25*cos(0.5*t);-2*cos(t)]; %Enter the expression

% Parse current states (X is the same size and "form" as X0)
% (i.e., in this sim, X = [e;r;thetahat])
e        = [X(1);X(2)];
r        = [X(3);X(4)];
thetaHat = [X(5);X(6);X(7);X(8);X(9)];

% Compute current x and xDot for convenience
q        = e + qd;
qDot     = r - a*e + qdDot;


% Compute cos(x2) and sin(x2) for convenience
c2       =cos(q(2));
s2       =sin(q(2));
c2d      =cos(qd(2));
s2d      =sin(qd(2));

% Compute current matrices for the dynamics
M        = [p1+2*p3*c2 p2+p3*c2;p2+p3*c2 p2];
Vm       = [-p3*s2*qDot(2) -p3*s2*(qDot(1)+qDot(2));p3*s2*qDot(1) 0];
fd       = [f1 0;0 f2];

eDot        = r - a*e;


mqqDotY=[qDot(1),qDot(2),2*c2*qDot(1)+c2*qDot(2),0,0;0,qDot(2)+qDot(1),c2*qDot(1),0,0];
mqdqdDotY=[qdDot(1),qdDot(2),2*c2d*qdDot(1)+c2d*qdDot(2),0,0;0,qdDot(2)+qdDot(1),c2d*qdDot(1),0,0];
fhDotd=ohmyd+beta*mqdqdDotY-beta*exp(-beta*t(length(t)))*Md0qdDotY;
fhDot=ohmy+beta*mqqDotY-beta*exp(-beta*t(length(t)))*M0qDotY;

fg=wy;
fgd=wyd;
Yf=fhDot+fg;
Yfd=fhDotd+fgd;
thetaTilde=theta-thetaHat;
Epsilon=uf-Yf*(thetaHat);


ohmyDot=-beta*ohmy-beta^2*mqqDotY;
ohmyMatDot=[ohmyDot(1,1:5)';ohmyDot(2,1:5)'];
ohmydDot=-beta*ohmyd-beta^2*mqdqdDotY;
ohmydMatDot=[ohmydDot(1,1:5)';ohmydDot(2,1:5)'];
mDotqDoty=[0,0,-2*s2*qDot(1)*qDot(2)-s2*qDot(2)*qDot(2) ,0 ,0;
           0,0,-s2*qDot(1)*qDot(2),                      0, 0];
mdDotqdDoty=[0,0,-2*s2d*qdDot(1)*qdDot(2)-s2d*qdDot(2)*qdDot(2) ,0 ,0;
             0,0,-s2d*qdDot(1)*qdDot(2),                         0, 0];
Ny=[0,0,-2*s2*qDot(2)*qDot(1)-s2*qDot(2)*qDot(2),qDot(1),0;
    0,0,s2*qDot(1)*qDot(1),0,qDot(2)];
Nyd=[0,0,-2*s2d*qdDot(1)*qdDot(2)-s2d*qdDot(2)*qdDot(2),qdDot(1),0;
     0,0,s2d*qdDot(1)*qdDot(1)                          ,0,qdDot(2)];
wyDot=-beta*wy+beta*(-mDotqDoty+Ny);
wyMatDot=[wyDot(1,1:5)';wyDot(2,1:5)'];
wydDot=-beta*wyd+beta*(-mdDotqdDoty+Nyd);
%disp(wydDot);
wydMatDot=[wydDot(1,1:5)';wydDot(2,1:5)'];
%disp(wydMatDot);
% Compute current regression matrix
y11      = a*eDot(1)-qdDotDot(1);%Enter the expression
y12      = a*eDot(2)-qdDotDot(2); %Enter the expression
y13      = 2*a*c2*eDot(1)+a*c2*eDot(2)-2*c2*qdDotDot(1)-c2*qdDotDot(2)-a*s2*qDot(2)*e(1)-a*s2*e(2)*qDot(1)-a*s2*e(2)*qDot(2)+s2*qDot(2)*qdDot(1)+s2*qDot(1)*qdDot(2)+s2*qDot(2)*qdDot(2); %Enter the expression
y14      = -qDot(1); %Enter the expression
y15      = 0; %Enter the expression
y21      = 0; %Enter the expression
y22      = a*eDot(1)+a*eDot(2)-qdDotDot(1)-qdDotDot(2); %Enter the expression
y23      = a*c2*eDot(1)-c2*qdDotDot(1)+a*s2*qDot(1)*e(1)-s2*qDot(1)*qdDot(1); %Enter the expression
y24      = 0; %Enter the expression
y25      = -qDot(2); %Enter the expression
Y        = [y11 y12 y13 y14 y15;y21 y22 y23 y24 y25];
 % Compute current regression matrix
y11      = qdDotDot(1);%Enter the expression
y12      = qdDotDot(2); %Enter the expression
y13      = 2*c2d*qdDotDot(1)+c2d*qdDotDot(2)-s2d*qdDot(2)*qdDot(1)-s2d*qdDot(2)*qdDot(2); %Enter the expression
y14      = qdDot(1); %Enter the expression
y15      = 0; %Enter the expression
y21      = 0; %Enter the expression
y22      = qdDotDot(1)+qdDotDot(2); %Enter the expression
y23      = c2d*qdDotDot(1)+s2d*qdDot(1)*qdDot(1); %Enter the expression
y24      = 0; %Enter the expression
y25      = qdDot(2); %Enter the expression
Yd        = [y11 y12 y13 y14 y15;y21 y22 y23 y24 y25];   
% Design controller

u        = -K*r-Y*thetaHat-e;%+Y*thetaHat;%-K*r; % -e Enter the expression
ufDot=-beta*uf+beta*u;

tilda=theta-thetaHat;
tau=[tau;t(length(t)),u'];
tilde=[tilde;t(length(t)),tilda'];
eps=[eps;t(length(t)),Epsilon'];
% Compute current closed-loop dynamics

rDot        = M\(-Vm*qDot-fd*qDot+u)-qdDotDot+a*eDot; %Enter the expression
thetaHatDot =gamma*eye(5)*transpose(Y)*r+gamma*eye(5)*transpose(Yf)*Epsilon;% eye(5)*transpose(Yd)*r;%Enter the expression


% Stacked dynamics vector (XDot is the same size and "form" as X)
XDot        = [eDot;rDot;thetaHatDot;ohmyMatDot;wyMatDot;ohmydMatDot;wydMatDot;ufDot];
end