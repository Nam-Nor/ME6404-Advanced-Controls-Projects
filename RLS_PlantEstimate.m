clear all
clc

T       = 0.04;                  % 40ms sample period
% syms a b c d e f g h z
% collect((a*z+b)/(z^2+c*z+d) + (e*z+f)/(z^2+g*z+h))

%% Use RLS to estimate discrete parameters of plant
u_1 = 0;u_2 = 0;u_3 = 0;u_4 = 0;
y_1 = 0;y_2 = 0;y_3 = 0;y_4 = 0;
m=0;n=0;o=0;p=0;q=0;r=0;s=0;t=0; %parameters to estimate

%needs to move 1.9199 radians from payload pickup to dropoff
time= 1.9199/0.35; % in s
N       = round(time/T);
Y       = zeros(N,1);
E       = Y;
TIME    = 0:T:T*(N-1);
THETA   = zeros(N,8);

P_1     = 20*eye(8);
Theta_1 = 0*[m n o p q r s t]';
lambda  = .8;

for k = 1:N

    %the real plant
    u    = 0.35; % slew velocity 
    Y(k) = Y(k)+u*T; % output is the angle of the payload
    y    = Y(k);
   
    %the parameter estimator
    Phi     = [u_1 u_2 u_3 u_4 -y_1 -y_2 -y_3 -y_4]';
    P       = (P_1 - (P_1*Phi * Phi' * P_1) / (1 + Phi' * P_1 * Phi)) / lambda;
    Theta   = Theta_1 + P * Phi * (y - (Phi' * Theta_1));
%     E(k)    = y - (Phi' * Theta_1);

    Theta_1 = Theta;
    P_1     = P;

    THETA(k,:) = Theta';
    
    y_4 = y_3;
    y_3 = y_2;
    y_2 = y_1;
    y_1 = y;
    u_4 = u_3;
    u_3 = u_2;
    u_2 = u_1;
    u_1 = u; 
        
end

%estimated parameters
m = Theta(1); 
n = Theta(2); 
o = Theta(3); 
p = Theta(4);
q = Theta(5); 
r = Theta(6); 
s = Theta(7); 
t = Theta(8);

%estimated system
num  = [m n o p];
den  = [1 q r s t];
sysd_e = tf(num, den, T);

%Controllable canonical realiation of SISO system
G=[0 1 0 0;0 0 1 0;0 0 0 1;-t -s -r -q];
H=[0 0 0 1]';
C=[p o n m];
D=0;

%Plant
plant=ss(G,H,C,D,T);
%sys_cont=d2c(sys_disc,'matched'); % convert to s domain

%% Model, estimated with zero damping estimate
g=9.81;
L1=0.821;  %suspension length (variable)
L2=0.7; % rigging length
mh=0.210; %hook mass
mp=0.075; %payload+tiny magnet mass
R=mp/mh;
beta=sqrt((1+R)^2*(1/L1+1/L2)^2 - 4*((1+R)/(L1*L2)));
w1=sqrt(g/2)*sqrt((1+R)*(1/L1+1/L2)+beta); %higher frequency
w2=sqrt(g/2)*sqrt((1+R)*(1/L1+1/L2)-beta);

sys1=tf(w1^2,[1 0 w1^2]);
sys2=tf(w2^2,[1 0 w2^2]);
model_cont=sys1+sys2;
model_disc=c2d(model_cont,T);
[num,den] = tfdata(model_disc);
num=cell2mat(num);
den=cell2mat(den);

Gm=[0 1 0 0;0 0 1 0;0 0 0 1;-den(5) -den(4) -den(3) -den(2)];
Hm=[0 0 0 1]';
Cm=[num(5) num(4) num(3) num(2)];
Dm=0;

model=ss(Gm,Hm,Cm,Dm,T);

%% *** Liapunov steps: pick Q and solve for P 
% Hint: use symbolic representation and the solve function
Q = eye(4);
syms p11 p22 p33 p44 p12 p13 p14 p23 p24 p34
P=[p11 p12 p13 p14;p12 p22 p23 p24;p13 p23 p33 p34;p14 p24 p34 p44];
Psolved=solve(-Q == transpose(Gm)*P+P*Gm);
Psolved=double(struct2array(Psolved));
% *** use P in the equation for M and solve for u, the control law
% P=[Psolved(1) Psolved(2);Psolved(2) Psolved(3)];
% p1=P(1,1);
% p2=P(1,2);
% p3=P(2,2);
% u=(v - wn_model^2*x1 - 2*zeta_model*wn_model*x2 +wn^2*x1 + 2*zeta*wn*x2)-M;

syms e1 e2 e3 e4 x1 x2 x3 x4 v u
X=[x1 x2 x3 x4]';

M=[e1 e2 e3 e4]*P*(Gm*X+Hm*v-G*X-H*u)



