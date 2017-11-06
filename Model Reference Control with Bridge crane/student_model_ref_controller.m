function u=student_model_ref_controller(xd1, xd2, v, x1, x2)

% declaring these variables global so we can see them here
global Ap  A  P;  %Bp  B ;
g=9.81;
L=1;   
zeta = 0.001;
wn = sqrt(g/L);
wn_model = wn;
zeta_model = 0.5; %Model A: 0.1, Model B: 0.5, Model C: 0.8

% errors
e1 = xd1 - x1;
e2 = xd2 - x2;
err=[e1 e2];
X=[x1; x2];
M=0;

% insert the control law here:
% u should be in m/s
% syms in
% eqn = err*P*(A*X-Ap*X-Bp*in+B*v)==M;
% law=solve(eqn,in);
% u=double(law);

p1=P(1,1);
p2=P(1,2);
p3=P(2,2);
u=(v - wn_model^2*x1 - 2*zeta_model*wn_model*x2 +wn^2*x1 + 2*zeta*wn*x2)-M;
return;