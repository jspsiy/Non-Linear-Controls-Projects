
 clear all
close all
%Set up parameters for sim
p1       = 3.473;
p2       = 0.196;
p3       = 0.242;
f1       = 5.3;
f2       = 1.1;
global tau sumtau2 tilde pos timeprev Yliprev Uiprev Y1hist sumu qDotprev sumY1i  qDotArr TIMER;
global c1Value qprev Y1iprev Y1i Mprev Vmprev MDotprev K KCL a gamma Y1imat count lambdamat tau2
global Y1eigenValMin lambda
lambda=1;
count=-1;
K        =[3,0;0,4];% [1,0;0,0.01];%9; %Enter a number
KCL      =0.00001;
a=0.3;
Y1eigenValMin=0;
gamma    = zeros(5);%75;
gamma(1,1)    = 30;
gamma(2,2)    = 1;%75;
gamma(3,3)    = 1;%75;
gamma(4,4)    = 100;%75;
gamma(5,5)    = 1;%75;
Y20=[0.6900,2.0000,5.2131,-2.2000,0;0,2.6900,-0.3621,0,0];
Uiprev=-K*[3;2]-Y20*[1;1;1;1;1]-[4,10];
TIMER=[0];
tau=[];
sumtau2=[0;0];
%Y1matadd
c1Value=[0;0;0;0;0];
tilde=[];
pos=0;
timeprev=0;
Y1hist=zeros(5);
qDotprev=[2.2;0];
Y1i=[];
qDotArr=[0;0];
qprev=[5;12];
Y1iprev=[0,0,0,2.2000,0;0,0,-2.5970,0,0];
%Mprev        = [p1+2*p3*1 p2+p3*1;p2+p3*1 p2];
Mprev        =[3.8814,0.4002;0.4002,0.1960];
MDotprev=[0 0;0 0];
Vmprev=[0,0.2857;-0.2857,0];
Y1imat=[];
tau2=[];
Ymat=[];
umat=[];
qDotmat=[];
Vmat=[];

% Stacked parameter vector
theta    = [p1;p2;p3;f1;f2];

% Simulation final time
tf   = 60;

% Initial condition vector (X0 must be same size and "form" as X and XDot below)
% (i.e., in this sim, X0 = [e0;r0;thetahat0])
X0   = [4;10;3;2;1;1;1;1;1;];

% Options for integration function
opts = odeset('OutputFcn',@dS,'RelTol',1e-3,'AbsTol',1e-3);

% Integrate (you can send the paramters theta to the dynamics as seen below)

[t,STATES] = ode45(@(t,X) twoLinkdynamics(t,X,theta),[0 tf],X0,opts);

% Set up desired trajectory data for plots (enter desired trajectory for your simulation)
qd = [cos(0.5*t) 2*cos(t)]';

% Parse integrated states (STATES is the same "form" as X0)
% (i.e., in this sim, STATES = [e r thetahat] over all time);
e  = STATES(:,1:2)';
r  = STATES(:,3:4)';
thetaHat = STATES(:,5:9)';
%input=STATES(:,10:11)';
% Compute x from e and xd for plotting purposes
q  = e + qd;

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
plot(lambdamat(:,1),lambdamat(:,2));

legend('minimum eigenvalue')
title("Eigenvalue plot")


function [XDot] = twoLinkdynamics(t,X,theta)
% initialize global variables for use in plot later.
global tau sumtau2 tilde pos timeprev Y1iprev Uiprev Y1hist sumu qDotprev sumY1i  qDotArr TIMER;
global c1Value qprev Y1i Mprev Vmprev MDotprev K a KCL gamma Y1imat count lambdamat tau2
global Y1eigenValMin
TIMER=[TIMER,t];

p1 = theta(1);
p2 = theta(2);
p3 = theta(3);
f1 = theta(4);
f2 = theta(5);

% Select gains for controller


% Desired trajectory and needed derivatives
qd       = [cos(0.5*t);2*cos(t)];
qdDot    = [-0.5*sin(0.5*t);-2*sin(t)]; %Enter the expression
qdDotDot = [-0.25*cos(0.5*t);-2*cos(t)]; %Enter the expression

% select gains (i.e., K, Kcl/Kicl, alpha, Gamma)
%missing%

% Parse current states (X is the same size and "form" as X0)
% (i.e., in this sim, X = [e;r;thetahat])
e        = [X(1);X(2)];
r        = [X(3);X(4)];
thetaHat = [X(5);X(6);X(7);X(8);X(9)];

% Compute current x and xDot for convenience
q        = e + qd;
qDot     = r - a*e + qdDot;
%qDotArr=[qDotArr,qDot];
eDot        = r - a*e;

% Compute cos(x2) and sin(x2) for convenience
c2       =cos(q(2));
s2       =sin(q(2));


% Compute current matrices for the dynamics
M        = [p1+2*p3*c2 p2+p3*c2;p2+p3*c2 p2];
MDot=[-2*p3*s2*qDot(2),-p3*s2*qDot(2);-p3*s2*qDot(2),0];
Vm       = [-p3*s2*qDot(2) -p3*s2*(qDot(1)+qDot(2));p3*s2*qDot(1) 0];
fd       = [f1 0;0 f2];


% compute current regression matrix (Y2 for gradient term/ non cl/icl)
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
 %Y2matadd=[y11 y12 y13 y14 y15,y21 y22 y23 y24 y25];
 Y2        = [y11 y12 y13 y14 y15;y21 y22 y23 y24 y25];

% set lambda for history stack min eigenvalue condition

if (pos==0)
    
end
% determine if the history stack meets the min eigenvalue condition
%if min eigenvalue of history stack less than lambda
    % compute qDotDot for Y1i
   % determine delta time.


    %if first loop (cannot determine dq or dt), initiallize qDotDot.
dt        = t(length(t)) - timeprev; 
dq        = qDot - qDotprev; % determine delta qDot.
qDotDot   = dq/dt; % determine qDotDot.
if dt==0
    qDotDot=[0;0];
end
    %determine Y1i (Y1 for this loop)
 y11      = qDotDot(1);%Enter the expression
 y12      = qDotDot(2); %Enter the expression
 y13      = c2*(2*qDotDot(1)+qDotDot(2))-s2*qDot(2)*(2*qDot(1)+ qDot(2)); %Enter the expression
 y14      = qDot(1); %Enter the expression
 y15      = 0; %Enter the expression
 y21      = 0; %Enter the expression
 y22      = qDotDot(1)+qDotDot(2); %Enter the expression
 y23      = c2*qDotDot(1)+s2*qDot(1)^2; %Enter the expression
 y24      = 0; %Enter the expression
 y25      = qDot(2); %Enter the expression
 Y1i        = [y11 y12 y13 y14 y15;y21 y22 y23 y24 y25]; 
% 
% 
%     % determine new history stack (summation)
% Y1hist    = Y1hist + Y1i'*Y1i;
%  
%     % optional: save the Y1i data in an array instead, then add/remove data
%     % from history stack to maximize min eigenvalue of Y1hist.
% 
%     % determine min eigenvalue of history stack
%     Y1eigenVal = eig(Y1hist);
% 
%     
%     % check min eigenvalue condition
%     Y1eigenValMin = min(Y1eigenVal);
%     
%     if (Y1eigenValMin > lambda)
%         
%        pos = 1;
%     else
%        pos = 0;
%     end
% end

% design controller (i.e., u)
u=-K*r-Y2*thetaHat-e;
tau=[tau;t(length(t)),u'];
%missing%
rDot        = M\(-Vm*qDot-fd*qDot+u)-qdDotDot+a*eDot; %Enter the expression

% determine which update law to use (i.e., thetaHatDot)
if pos == 1
    %(after eigenvalues condition met for Y1hist)
    %if Y1eigenValMin>max(lamb
    %Y1imat        = [Y1imat;y11 y12 y13 y14 y15 y21 y22 y23 y24 y25];
    %u2=M*qDotDot+Vm*qDot+F;
    %tau2=[tau2;u'];    
    % for CL:
    % compute new cl_sumY (function of: Yi,ui,thetaHati),
   %thetaHatDot =(gamma)*transpose(Y2)*r+KCL*gamma*sumY1i'*(sumtau2-sumY1i*thetaHat);                                     %c1Value=c1Value+Y1i'*(u-Y1i*thetaHat);
    if Y1eigenValMin>max(lambdamat)
      Y1imat        = [Y1imat;Y1i(1,:),Y1i(2,:)];
     %u2=M*qDotDot+Vm*qDot+F;
     tau2=[tau2;u'];
    cl_sumY(Y1imat,tau2,thetaHat); 
    end
   thetaHatDot =(gamma)*transpose(Y2)*r+KCL*gamma*c1Value;%1/length(Y1imat(:,1))*KCL*gamma*c1Value;
    %if Y1eigenValMin>max(lambdamat)
    %Y1imat        = [Y1imat;y11 y12 y13 y14 y15 y21 y22 y23 y24 y25];
    %tau2=[tau2;u'];
%    cl_sumY(Y1imat,tau2,thetaHat);
    %end
    % assuming new data in history stack
    % for ICL:
    
    pp=0;
    % compute new icl_sumY(script) (function of: Yi(script),ui,x(ti),x(ti-dt),thetaHati),
    % assuming new data in history stack

    % compute thetaHatDot
    %missing%
else
    Y1imat        = [Y1imat;Y1i(1,:),Y1i(2,:)];
    %u2=M*qDotDot+Vm*qDot+F;
    tau2=[tau2;u'];
    % (before eigenvalues condition met for Y1hist)

    thetaHatDot =(gamma)*transpose(Y2)*r;    % for CL:
    % compute cl_sumY (function of: Yi,ui, thetaHat) for use later

     cl_sumY(Y1imat,tau2,thetaHat);
   %c1Value=c1Value+Y1i'*(u-Y1i*thetaHat);
%c1Value=cl_sumY(Yi,ui, thetaHat);
    % for ICL:
    % compute icl_sumY(script) (function of: Yi(script),ui,x(ti),x(ti-deltat),thetaHati) for use later
    %missing%
    
end 

% compute current closed-loop errors for integration(i.e., eDot, rDot)
%missing%

% update "previous" variables for next loop
%(i.e., qDotprev, timeprev, Y1iprev, Uiprev)
%missing%
qDotprev= qDot;
qprev=q;

Y1iprev=Y1i;
Uiprev=u;
Mprev=M;
MDotprev=MDot;
Vmprev=Vm;
tilda=theta-thetaHat;
tilde=[tilde;t,tilda'];
lambdamat=[lambdamat;t,Y1eigenValMin];
% Stacked dynamics vector (XDot is the same size and "form" as X)
XDot        = [eDot;rDot;thetaHatDot];

end

function status= dS(t,y,flag)
global a timeprev qDotprev Y1hist lambda pos lambdamat Y1eigenValMin Y1i count
count=count+1;
%if (length(t)>0 & (sum(t==60)~=1 ))
if (count>0 & length(t)>0)
qd = [cos(0.5*t(length(t)));2*cos(t(length(t)))]';
qdDot=[-0.5*sin(0.5*t(length(t)));-2*sin(t(length(t)))];

e        = [y(1);y(2)];
r        = [y(3);y(4)];   
q        = e + qd;
c2=cos(q(2));
s2=sin(q(2));
qDot     = r - a*e + qdDot;
dt        = t(length(t)) - timeprev; 
dq        = qDot - qDotprev; % determine delta qDot.
qDotDot   = dq/dt; % determine qDotDot.
if dt==0
    qDotDot=[0;0];
end
timeprev=t(length(t));    
y11      = qDotDot(1);%Enter the expression
y12      = qDotDot(2); %Enter the expression
y13      = c2*(2*qDotDot(1)+qDotDot(2))-s2*qDot(2)*(2*qDot(1)+ qDot(2)); %Enter the expression
y14      = qDot(1); %Enter the expression
y15      = 0; %Enter the expression
y21      = 0; %Enter the expression
y22      = qDotDot(1)+qDotDot(2); %Enter the expression
y23      = c2*qDotDot(1)+s2*qDot(1)^2; %Enter the expression
y24      = 0; %Enter the expression
y25      = qDot(2); %Enter the expression
Y1i        = [y11 y12 y13 y14 y15;y21 y22 y23 y24 y25]; 


    % determine new history stack (summation)
Y1hist    = Y1hist + Y1i'*Y1i;
 
    % optional: save the Y1i data in an array instead, then add/remove data
    % from history stack to maximize min eigenvalue of Y1hist.

    % determine min eigenvalue of history stack
    Y1eigenVal = eig(Y1hist);

    
    % check min eigenvalue condition
    Y1eigenValMin = min(Y1eigenVal);
    
    if (Y1eigenValMin > lambda)
        
       pos = 1;
     else
       pos = 0;
%      %  lambdamat=[lambdamat;t,Y1eigenValMin];
    end
end
status=0;    
end

function cl_sumY(Y1imat,tau2,thetaHat)
global c1Value
   c1Value=0;
   
   for i=1:length(tau2(:,1))
       Yi=[Y1imat(i,1),Y1imat(i,2),Y1imat(i,3),Y1imat(i,4),Y1imat(i,5);
           Y1imat(i,6),Y1imat(i,7),Y1imat(i,8),Y1imat(i,9),Y1imat(i,10)];
       ui=[tau2(i,1);tau2(i,2)];
       c1Value=c1Value+Yi'*(ui-Yi*thetaHat);
   end
    end
