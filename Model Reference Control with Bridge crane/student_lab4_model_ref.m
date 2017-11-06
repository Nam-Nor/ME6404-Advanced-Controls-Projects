%=========================================================
%      ME6404: ADVANCED CONTROLS AND IMPLEMENTATION
%                       
%             LAB 4:  MODEL REFERENCE CONTROL
%==========================================================
clear;clc;

%% Part A: Model Reference control on the Trolley Velocity
%  PARAMETERS THAT DEFINE THE MODEL  (CHANGE THESE)
tau_model=1;    %first order system time constant
zeta_model=1;   %second order system damping ratio
w_model=1;     %second order system natural frequency
 
% FORM THE CONTINUOUS-TIME STATE SPACE MATRICES 
% Hint: use the ss command
 

%  INSERT CODE HERE TO RUN A TIME SIMULATION OF YOUR
%  FIRST- AND SECOND-ORDER MODELS TO PULSE INPUTS
%  Hint: use the lsim command
 


%% Part B: Model Reference Control on Payload Oscillation

% **** make a pulse velocity command signal, v
% velocity step size of the portable bridge crane
v_step = 0.243;      % 0.243m/s is 100% speed on the PBC
step1 = v_step;
step2 = -v_step;

% times of the velocity input pulse signal, v
t_on = 1;
t_off = 5;
t_sim = 10;         % simulation time

% declaring these variables global
global Ap Bp Cp Dp A B C D P;

% *** Plant (Payload)
g=9.81;
L=1;                      %suspension length
zeta = 0.001;
wn = sqrt(g/L);
Ap = [0 1; -wn^2 -2*zeta*wn];
Bp = [0; 1];
Cp = eye(2);
Dp = [0;0];

% *** Model
wn_model = wn;
zeta_model = 0.5; %Model A: 0.1, Model B: 0.5, Model C: 0.8
A = [0 1; -wn_model^2 -2*zeta_model*wn_model];
B = [0; 1];
C = eye(2);
D = [0;0];

% *** Liapunov steps: pick Q and solve for P 
% Hint: use symbolic representation and the solve function
Q = eye(2);
syms p1 p2 p3
P=[p1 p2;p2 p3];
Psolved=solve(-Q == transpose(A)*P+P*A);
Psolved=double(struct2array(Psolved));
% *** use P in the equation for M and solve for u, the control law
P=[Psolved(1) Psolved(2);Psolved(2) Psolved(3)];
%% run the continuous simulink model
sim('student_model_ref');
