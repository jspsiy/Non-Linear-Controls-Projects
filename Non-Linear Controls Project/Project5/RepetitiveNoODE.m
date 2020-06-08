 clear all
close all
%Set up parameters for sim
p1       = 3.473;
p2       = 0.196;
p3       = 0.242;
f1       = 5.3;
f2       = 1.1;
fs1=1.2;
fs2=0.4;
fs=[fs1,0;0,fs2];
global tau tilde;
tau=[];
tilde=[];
wHatprev=0;
wHatmat=[];
wmat=[];
timestack=[];
% Stacked parameter vector
theta    = [p1;p2;p3;f1;f2];
zeta=[fs1;fs2];
% Simulation final time
tf   = 60;

% Initial condition vector (X0 must be same size and "form" as X and XDot below)
% (i.e., in this sim, X0 = [e0;r0;thetahat0])
X0   = [4;10;3;2;0;0;];
Xstates=[];
Xstates=[Xstates;X0'];
% Options for integration function
opts = odeset('RelTol',1e-3,'AbsTol',1e-3);

% Integrate (you can send the paramters theta to the dynamics as seen below)

% Select gains for controller
K        = [5,0;0,2];%9; %Enter a number
KL1      =3;
KL2       =3;
a        = 0.8;%0.15; %Enter a number
gamma    = zeros(2);%75;
gamma(1,1)    = 50;%75;
gamma(2,2)    = 50;%75;;
beta=3;
count=1;
dt=0.001;

X=X0;
for t=0:dt:tf
    
w1=0.5;   w2=1;
T1=round(2*pi/w1/dt)+1;
T2=round(2*pi/w2/dt)+1;
timestack=[timestack;t];
% Desired trajectory and needed derivatives
qd       = [cos(w1*t);2*cos(w2*t)];
qdDot    = [-0.5*sin(w1*t);-2*sin(w2*t)]; %Enter the expression
qdDotDot = [-0.25*cos(w1*t);-2*cos(w2*t)]; %Enter the expression

% Parse current states (X is the same size and "form" as X0)
% (i.e., in this sim, X = [e;r;thetahat])
e        = [X(1);X(2)];
r        = [X(3);X(4)];
zetaHat = [X(5);X(6)];

% Compute current x and xDot for convenience
q        = qd - e;
qDot     = -r + a*e + qdDot;


% Compute cos(x2) and sin(x2) for convenience
c2       =cos(q(2));
s2       =sin(q(2));
c2d      =cos(qd(2));
s2d      =sin(qd(2));

% Compute current matrices for the dynamics
M        = [p1+2*p3*c2 p2+p3*c2;p2+p3*c2 p2];
Vm       = [-p3*s2*qDot(2) -p3*s2*(qDot(1)+qDot(2));p3*s2*qDot(1) 0];
fd       = [f1 0;0 f2];
Mqd       =[p1+2*p3*c2d p2+p3*c2d;p2+p3*c2d p2];
Vmqd= [-p3*s2d*qdDot(2) -p3*s2d*(qdDot(1)+qdDot(2));p3*s2d*qdDot(1) 0];

w=Mqd*qdDotDot++Vmqd*qdDot+fd*qdDot;
eDot        = r - a*e;
% Compute current regression matrix
y11      = sign(qDot(1)); %Enter the expression
y22      = sign(qDot(2)); %Enter the expression
Y        = [y11,0;0,y22];

% Design controller


if count<=T1
    wHat1=KL1*r(1);
else
wHat1=saturation(beta,wHatmat(:,2),count,T1)+KL1*r(1);
end


if count<=T2
    wHat2=KL2*r(2);
else
wHat2=saturation(beta,wHatmat(:,3),count,T2)+KL2*r(2);
end

wHat=[wHat1;wHat2];



u        = K*r+Y*zetaHat+e+wHat;%+wHat;%+Y*thetaHat;%-K*r; % -e Enter the expression

zilda=zeta-zetaHat;
tau=[tau;t(length(t)),u'];
tilde=[tilde;t(length(t)),zilda'];

%plot(t(length(t)),u);
% Compute current closed-loop dynamics

rDot        =qdDotDot-M\(-Vm*qDot-fd*qDot-fs*sign(qDot)+u)+a*eDot; %Enter the expression
zetaHatDot =(gamma)*transpose(Y)*r;% eye(5)*transpose(Yd)*r;%Enter the expression

% Stacked dynamics vector (XDot is the same size and "form" as X)
XDot        = [eDot;rDot;zetaHatDot];
X=X+XDot*dt;
Xstates=[Xstates;X'];

wHatmat=[wHatmat;t,wHat(1),wHat(2)];
wmat=[wmat;t,w(1),w(2)];
wHatprev=wHat;
count=count+1;
end


Xstates=Xstates(1:(size(Xstates(:,1))-1),:);
% Set up desired trajectory data for plots (enter desired trajectory for your simulation)
qd = [cos(0.5*timestack) 2*cos(timestack)];

% Parse integrated states (STATES is the same "form" as X0)
% (i.e., in this sim, STATES = [e r thetahat] over all time);
e  = Xstates(:,1:2);
%e(end)=[];

r  = Xstates(:,3:4);

thetaHat =Xstates(:,5:6)';
%input=STATES(:,10:11)';
% Compute x from e and xd for plotting purposes
q  = qd - e;

% Plot the actual vs desired trajectories
subplot(2,2,1)
plot(timestack,qd(:,1),'-','LineWidth',2)
hold on
ax = gca;
ax.ColorOrderIndex = 1;
plot(timestack,q(:,1),':','LineWidth',2)
plot(timestack,qd(:,2),'-','LineWidth',2)
plot(timestack,q(:,2),':','LineWidth',2)

title('q vs qd');
hold off

% Plot the filtered tracking error
subplot(2,2,2)
hold on
plot(timestack,r(:,1),'--','LineWidth',2)
plot(timestack,r(:,2),'--','LineWidth',2)
plot(timestack,e(:,1),':','LineWidth',2)
plot(timestack,e(:,2),':','LineWidth',2)
legend("r1","r2","e1","e2");
title('Error')
hold off

% Plot the adaptive estimates vs actual parameters
subplot(2,2,3)
plot(timestack,repmat(zeta,1,length(timestack)),'-','LineWidth',2)
hold on
ax = gca;
ax.ColorOrderIndex = 1;
plot(timestack,thetaHat(1,:),':','LineWidth',2)
plot(timestack,thetaHat(2,:),':','LineWidth',2)
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

legend('tilde1','tilde2')
title("thetatilde")


figure(3)
plot(timestack,wHatmat(:,2),':','LineWidth',2)
hold on
plot(timestack,wHatmat(:,3),':','LineWidth',2)
title('wHat');
hold off


figure(4)
plot(timestack,wmat(:,2)-wHatmat(:,2),':','LineWidth',2)
hold on
plot(timestack,wmat(:,3)-wHatmat(:,3),':','LineWidth',2)
title('wTilde');
hold off


function [value] = saturation(beta,wHatmat,count,T)


wHatprev=wHatmat(count-T);
if wHatprev>=beta
    value=beta;
elseif wHatprev<=-beta
    value=-beta;
else
    value=wHatprev;
end
end