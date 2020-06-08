%function twoLinkRobotAdaptive
close all
clear all
%Set up parameters for sim
p1       = 3.473;
p2       = 0.196;
p3       = 0.242;
f1       = 5.3;
f2       = 1.1;
global tau tilde eps e20 thetahat0 thetahatplot;
tau=[];
tilde=[];
eps=[];
% Stacked parameter vector
theta    = [p1;p2;p3;f1;f2];

% Simulation final time
tf   = 60;

% Initial condition vector (X0 must be same size and "form" as X and XDot below)
% (i.e., in this sim, X0 = [e0;r0;thetahat0])
thetahat0=[1;1;1;1;1];
X0   = [4;10;%e1 @2
        3;2;%e2  @4
        %1;1;1;1;1;%thetahat @9
        0;0;0;0;0;%thetahat2 called zetahat @9
        0;0;0;0;0; %Yfd states @14
        0;0;0;0;0; %Yfd states @19

        0;0;  % ufd states @31
        0;0; %mew states @33
        0;0]; % mew2 states @35
% Options for integration function
opts = odeset('RelTol',1e-3,'AbsTol',1e-3);

qd0 = [cos(0.5*0);2*cos(0)];
qdDot0 = [-0.5*sin(0.5*0);-2*sin(0)];
e10=[X0(1);X0(2)];
e20=[X0(3);X0(4)];
q0=qd0-e10;
a=0.2;
%qDot0=qdDot0-a1*e10-e20;

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
zetaHat = STATES(:,5:9)';

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
%plot(t,r,'--','LineWidth',2)
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
plot(thetahatplot(:,1),thetahatplot(:,2))
plot(thetahatplot(:,1),thetahatplot(:,3))
plot(thetahatplot(:,1),thetahatplot(:,4))
plot(thetahatplot(:,1),thetahatplot(:,5))
plot(thetahatplot(:,1),thetahatplot(:,6))
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
global tau tilde eps e20 thetahat0 thetahatplot;
% Parse parameter vector
p1 = theta(1);
p2 = theta(2);
p3 = theta(3);
f1 = theta(4);
f2 = theta(5);
e1       = [X(1);X(2)];
e2       = [X(3);X(4)];
zetaHat = [X(5);X(6);X(7);X(8);X(9)];
Yfd=[X(10:14)';X(15:19)'];
%omega=[X(20:24)';X(25:29)'];%omega parameter
%wyd=[X(20:24)';X(25:29)'];
uf=[X(20);X(21)];
mew=[X(22);X(23)];
mew2=[X(24);X(25)];
% Select gains for controller
K1        =50; %[15,0;0,1.5];%9; %Enter a number
K2       =15;
K=[K1 0; 0 K2];
a1        =1.5;%0.15; %Enter a number
a2        =1.2;
gamma    = zeros(5);%75;
gamma(1,1)    = 500;%75;
gamma(2,2)    = 10;%75;
gamma(3,3)    = 3;%75;
gamma(4,4)    =300;%75;
gamma(5,5)    = 10;%75;
% gamma    = zeros(5);%75;
% gamma(1,1)    =30;%75;
% gamma(2,2)    =;%75;
% gamma(3,3)    =0.3;%75;
% gamma(4,4)    =20;%75;
% gamma(5,5)    =1;%75;
%gamma=5;
beta1=0.5;
%beta1=3;
beta2=0.9;
% Desired trajectory and needed derivatives
qd       = [cos(0.5*t);2*cos(t)];
qdD    = [-0.5*sin(0.5*t);-2*sin(t)]; %Enter the expression
qdDot=qdD;
qdDD = [-0.25*cos(0.5*t);-2*cos(t)]; %Enter the expression
qdDotDot=qdDD;
qdDDD=[0.125*sin(0.5*t);2*sin(t)];
qdDotDotDot=qdDDD;
qdDDDD = [0.0625*cos(0.5*t);2*cos(t)];
qdDotDotDotDot=qdDDDD;
% Parse current states (X is the same size and "form" as X0)
% (i.e., in this sim, X = [e;r;thetahat])


% Compute current x and xDot for convenience
q        = qd-e1;
e1Dot    =e2-a1*e1;
%e2Dot    = r - a2*e2 ;
qDot=qdD+a1*e1-e2;
%qDotDot=qdDotDot-r+a2*e2+a1*e1Dot;
%e1DotDot=r-a1*e1Dot-a2*e2;
% Compute cos(x2) and sin(x2) for convenience
c2       =cos(q(2));
s2       =sin(q(2));
c2d      =cos(qd(2));
s2d      =sin(qd(2));
cd2=c2d;
sd2=s2d;

% Compute current matrices for the dynamics
M        = [p1+2*p3*c2, p2+p3*c2;p2+p3*c2 ,p2];
%mDot=[-2*p3*s2*qDot(2),-p3*s2*qDot(2);-p3*s2*qDot(2),0];
Vm       = [-p3*s2*qDot(2) ,-p3*s2*(qDot(1)+qDot(2));p3*s2*qDot(1), 0];
%VmDot=[-p3*c2*qDot(2)^2-p3*s2*qDotDot(2),-p3*c2*qDot(2)*(qDot(1)+qDot(2))-p3*s2*(qDotDot(1)+qDotDot(2));p3*c2*qDot(2)*qDot(1)+p3*s2*qDotDot(1),0];
fd       = [f1 0;0 f2];
tauD=[0.5*cos(0.5*t);sin(t)];
%tauDDot=[-0.25*sin(0.5*t);cos(t)];



%mqqDotY=[qDot(1),qDot(2),2*c2*qDot(1)+c2*qDot(2),0,0;0,qDot(2)+qDot(1),c2*qDot(1),0,0];
%mqdqdDotY=[qdD(1),qdD(2),2*c2d*qdD(1)+c2d*qdD(2),0,0;0,qdD(2)+qdD(1),c2d*qdD(1),0,0];
%fhDotd=ohmyd+beta*mqdqdDotY-beta*exp(-beta*t(length(t)))*Md0qdDotY;
%fhDot=ohmy+beta*mqqDotY-beta*exp(-beta*t(length(t)))*M0qDotY;

%fg=wy;
%fgd=wyd;
%Yf=fhDot+fg;
%Yfd=fhDotd+fgd;


%ohmyDot=-beta*ohmy-beta^2*mqqDotY;
%ohmyMatDot=[ohmyDot(1,1:5)';ohmyDot(2,1:5)'];
%ohmydDot=-beta*ohmyd-beta^2*mqdqdDotY;
%ohmydMatDot=[ohmydDot(1,1:5)';ohmydDot(2,1:5)'];
%mDotqDoty=[0,0,-2*s2*qDot(1)*qDot(2)-s2*qDot(2)*qDot(2) ,0 ,0;
%           0,0,-s2*qDot(1)*qDot(2),                      0, 0];
%mdDotqdDoty=[0,0,-2*s2d*qdD(1)*qdD(2)-s2d*qdD(2)*qdD(2) ,0 ,0;
%             0,0,-s2d*qdD(1)*qdD(2),                         0, 0];
%Ny=[0,0,-2*s2*qDot(2)*qDot(1)-s2*qDot(2)*qDot(2),qDot(1),0;
%    0,0,s2*qDot(1)*qDot(1),0,qDot(2)];
%Nyd=[0,0,-2*s2d*qdD(2)*qdD(1)-s2d*qdD(2)*qdD(2),qdD(1),0;
%     0,0,s2d*qdD(1)*qdD(1)                          ,0,qdD(2)];
%wyDot=-beta*wy+beta*(-mDotqDoty+Ny);
%wyMatDot=[wyDot(1,1:5)';wyDot(2,1:5)'];
%wydDot=-beta*wyd+beta*(-mdDotqdDoty+Nyd);
%disp(wydDot);
%wydMatDot=[wydDot(1,1:5)';wydDot(2,1:5)'];
%disp(wydMatDot);
% Compute current regression matrix
%  y11      = a1*e1Dot(1)-qdDotDot(1);%Enter the expression
%  y12      = a1*e1Dot(2)-qdDotDot(2); %Enter the expression
%  y13      = 2*a1*c2*e1Dot(1)+a1*c2*e1Dot(2)-2*c2*qdDotDot(1)-c2*qdDotDot(2)-a1*s2*qDot(2)*e1(1)-a1*s2*e1(2)*qDot(1)-a1*s2*e1(2)*qDot(2)+s2*qDot(2)*qdDot(1)+s2*qDot(1)*qdDot(2)+s2*qDot(2)*qdDot(2); %Enter the expression
%  y14      = -qDot(1); %Enter the expression
%  y15      = 0; %Enter the expression
%  y21      = 0; %Enter the expression
%  y22      = a1*e1Dot(1)+a1*e1Dot(2)-qdDotDot(1)-qdDotDot(2); %Enter the expression
%  y23      = a1*c2*e1Dot(1)-c2*qdDotDot(1)+a1*s2*qDot(1)*e1(1)-s2*qDot(1)*qdDot(1); %Enter the expression
%  y24      = 0; %Enter the expression
%  y25      = -qDot(2); %Enter the expression
%  Y        = [y11 y12 y13 y14 y15;y21 y22 y23 y24 y25];
 % Compute current regression matrix
y11      = qdDD(1);%Enter the expression
y12      = qdDD(2); %Enter the expression
y13      = c2d*(2*qdDD(1)+qdDD(2))-s2d*qdD(2)*(2*qdD(1)+ qdD(2)); %Enter the expression
y14      = qdD(1); %Enter the expression
y15      = 0; %Enter the expression
y21      = 0; %Enter the expression
y22      = qdDD(1)+qdDD(2); %Enter the expression
y23      = c2d*qdDD(1)+s2d*qdD(1)^2; %Enter the expression
y24      = 0; %Enter the expression
y25      = qdD(2); %Enter the expression
Yd        = [y11 y12 y13 y14 y15;y21 y22 y23 y24 y25];   

Yd_0 =[-0.2500,-2.0000,1.0404,0,0,;0,-2.2500,0.1040,0,0];

%Compute Regression matrix for YdDot
y11      = qdDDD(1);%Enter the expression
y12      = qdDDD(2); %Enter the expression
%y13      = -2*s2d*qdD(2)*qdDD(1)+2*c2d*qdDDD(1)-s2d*qdD(2)*qdDD(2)+c2d*qdDDD(2)-c2d*qdD(2)^2*qdD(1)-s2d*qdDD(2)*qdD(1)-c2d*qdD(2)^3-2*s2d*qdD(2)*qdDD(2);%-s2d*qdDot(2)*qdDot(2); %Enter the expression
y13      =-2*s2d*qdD(2)*(2*qdDD(1)+qdDD(2))+c2d*(2*qdDDD(1)+qdDDD(2))-(c2d*qdD(2)^2+s2d*qdDD(2))*(2*qdD(1)+qdD(2));
%ydd13=cd2*(2*qdDotDotDot(1)+qdDotDotDot(2))-(2*qdDot(1)+qdDot(2))*(cd2*qdDot(2)*qdDot(2)+sd2*qdDotDot(2))-2*sd2*qdDot(2)*(2*qdDotDot(1)+qdDotDot(2));
y14      = qdDD(1); %Enter the expression
y15      = 0; %Enter the expression
y21      = 0; %Enter the expression
y22      = qdDDD(1)+qdDDD(2); %Enter the expression
y23      = -s2d*qdD(2)*qdDD(1)+c2d*qdDDD(2)            +c2d*qdD(2)*qdD(1)^2+2*s2d*qdDD(1)*qdD(1);%+s2d*qdDot(1)*qdDot(1); %Enter the expression
%ydd23    = -sd2*qdDot(2)*qdDotDot(1)+cd2*qdDotDotDot(2)+c2d*qdDot(1)^2*qdDot(2)-qdDot(2)*sd2*qdDotDot(2);
y24      = 0; %Enter the expression
y25      = qdDD(2); %Enter the expression
YdDot       = [y11 y12 y13 y14 y15;y21 y22 y23 y24 y25];   

YdDot_0=[0,0,0,-.25,0;0,0,0,0,-2];

%Compute Regression matrix for YdDotDot
y11      = qdDDDD(1);%Enter the expression
y12      = qdDDDD(2); %Enter the expression
%y13      = -2*s2d*qdD(2)*qdDD(1)+2*c2d*qdDDD(1)-s2d*qdD(2)*qdDD(2)+c2d*qdDDD(2)-c2d*qdD(2)^2*qdD(1)-s2d*qdDD(2)*qdD(1)-c2d*qdD(2)^3-2*s2d*qdD(2)*qdDD(2);%-s2d*qdDot(2)*qdDot(2); %Enter the expression
y13      =-2*s2d*qdDD(2)*(2*qdDD(1)+qdDD(2))-2*c2d*qdD(2)^2*(2*qdDD(1)+qdDD(2))-3*s2d*qdD(2)*(2*qdDDD(1)+qdDDD(2))+c2d*(2*qdDDDD(1)+qdDDDD(2))+2*s2d*qdD(2)^3*qdD(1)-4*c2d*qdD(2)*qdDD(2)*qdD(1)-2*c2d*qdD(2)^2*qdDD(1)-2*c2d*qdD(2)*qdD(1)*qdDD(2)-2*s2d*qdDDD(2)*qdD(1)-2*s2d*qdDD(2)*qdDD(1)+s2d*qdD(2)^4-3*c2d*qdD(2)^2*qdDD(2)-c2d*qdD(2)^2*qdDD(2)-s2d*qdDDD(2)*qdD(2)-s2d*qdDD(2)^2;           
y14      = qdDDD(1); %Enter the expression
y15      = 0; %Enter the expression
y21      = 0; %Enter the expression
y22      = qdDDDD(1)+qdDDDD(2); %Enter the expression
%y23      = -c2d*qdD(2)^2*qdDD(1)-s2d*qdDD(2)*qdDD(1)-2*s2d*qdD(2)*qdDDD(1)-s2d*qdD(2)*qdDD(1)+c2d*qdDDDD(1)-s2d*qdD(1)^2*qdD(2)^2+2*c2d*qdD(1)*qdDD(1)*qdD(2)+c2d*qdD(1)^2*qdDD(2)+2*c2d*qdD(2)*qdD(1)*qdDD(1)+2*s2d*qdD(1)*qdDDD(1);%+s2d*qdDot(1)*qdDot(1); %Enter the expression
y23       = -c2d*qdD(2)^2*qdDD(1)-s2d*qdDD(2)*qdDD(1)-s2d*qdD(2)*qdDDD(1)-s2d*qdD(2)*qdDD(1)+c2d*qdDDDD(1)-s2d*qdD(1)^2*qdD(2)^2+2*c2d*qdD(1)*qdDD(1)*qdD(2)+c2d*qdD(1)^2*qdDD(2)+2*c2d*qdD(2)*qdD(1)*qdDD(1)+2*s2d*qdDD(1)^2+2*s2d*qdD(1)*qdDDD(1);
y24      = 0; %Enter the expression
y25      = qdDDD(2); %Enter the expression
YdDotDot       = [y11 y12 y13 y14 y15;y21 y22 y23 y24 y25];   



% Design controller

thetaHat=thetahat0+gamma*(YdDot)'*e2-gamma*(YdDot_0)'*e20+zetaHat;

u        = Yd*thetaHat+(K+1)*(e2-e20)+mew;
%u        = mew;%((K+1)*e2 - (K+1)*e20 + mew);%+Y*thetaHat;%-K*r; % -e Enter the expr..ession
e2Dot=qdDD-M\(-Vm*qDot-fd*qDot-tauD+u)+a1*e1Dot;
%r=e2Dot+a2*e2;
muDot=(K+eye(2))*a2*e2 + beta1*sign(e2);%(K+1)*a2*e2 + beta*sign(e2);%+Y*thetaHat;%+Y*thetaHat;%-K*r; % -e Enter the expr..ession
ufHat=Yfd*thetaHat+mew2;
Epsilon=uf-ufHat;
omegaDot=zeros(2,5);
%omegaDot=-omega.*beta2-Yd.*(beta2^2);
%YfdDot=omega+Yd.*beta2-Yd_0.*(beta2*exp(-beta2*t));
YfdDot= -beta1*Yfd+beta1*Yd;
%thetaHatDot =gamma*eye(5)*transpose(YdDot)*(e2Dot+a2*e2)+gamma*eye(5)*transpose(YfdDot)*Epsilon;
zetaHatDot=-gamma*(e2'*YdDotDot)'+a2*gamma*YdDot'*e2+gamma*YfdDot'*Epsilon;
mew2Dot=K2*Epsilon+beta1*sign(Epsilon);
ufdDot=-beta2*uf+beta2*u;
%e2Dot=qdDD-M\(-Vm*qDot-fd*qDot-tauD+u)+a1*e1Dot;%+a2*e2;

%tauDot=YdDot*thetaHat+Yd*thetaHatDot+muDot;
%S=M*qdDotDot+Vm*qDot+fd*qDot+M*(a1*e1Dot+a2*e2)-Yd*theta;
%SDot=mDot*qdDotDot+M*qdDotDotDot+VmDot*qDot+Vm*qDotDot+fd*qDotDot+mDot*(a1*e1Dot+a2*e2)+M*(a1*e1DotDot+a2*e2Dot)-YdDot*theta+Yd*thetaHatDot;

tilda=theta-thetaHat;
tau=[tau;t(length(t)),u'];
tilde=[tilde;t(length(t)),tilda'];
eps=[eps;t(length(t)),Epsilon'];
thetahatplot=[thetahatplot;t(length(t)),thetaHat'];



% Stacked dynamics vector (XDot is the same size and "form" as X)
XDot        = [e1Dot;e2Dot;zetaHatDot;YfdDot(1,:)';YfdDot(2,:)';ufdDot;muDot;mew2Dot];
end