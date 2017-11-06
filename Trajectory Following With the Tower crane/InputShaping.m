clear;
clc;

% [65,567,1720] for start location
% 821 mm max hoist position for obstacle clearance
% [307,868,1580] for payload pickup; 1561 min hoist pos for magnet engage
% 879 mm max hoist poisition for payload pickup
% max trolley position for obstacle avoidance
% [197,811,835] for payload dropoff location
% 988 mm hoist position for payload dropoff

g=9.81;
L1=0.821;  %suspension length (variable)
L2=0.7; % rigging length
zeta = 0.001; %Damping ratio estimate
mh=0.210; %hook mass
mp=0.075; %payload+tiny magnet mass

R=mp/mh;
beta=sqrt((1+R)^2*(1/L1+1/L2)^2 - 4*((1+R)/(L1*L2)));
w1=sqrt(g/2)*sqrt((1+R)*(1/L1+1/L2)+beta); %higher frequency
w2=sqrt(g/2)*sqrt((1+R)*(1/L1+1/L2)-beta);

t1=2*pi/w1;
t2=2*pi/w2;
Vt=0.20; %tolerable vibration
side=(1+Vt)/4;
middle=(1-Vt)/2;

%First mode ZV shaper
ZV1=[0.5, 0.5;0, (t1)/2];
%First mode EI shaper
EI1=[side,middle,side; 0,(t1)/2, t1];
%Second mode ZV shaper
ZV2=[0.5, 0.5;0, (t2)/2];
%Second mode EI shaper
EI2=[side,middle,side; 0,(t2)/2, t2];
%Dual mode ZV shaper
ZVdual=[0.25 0.25 0.25 0.25; 0 t1/2 t2/2 (t1+t2)/2]
%Dual mode EI shaper
shaper1=EI1;
shaper2=EI2;
s1Len = length(shaper1);
s2Len = length(shaper2);
shaperConv = zeros(2,s1Len*s2Len);
for i = 1:s1Len
    for j = 1:s2Len
        v = shaper1(1,i)*shaper2(1,j);
        t = shaper1(2,i) + shaper2(2,j);
        
        ind = s1Len*(i-1)+j;
        shaperConv(1,ind) = v;
        shaperConv(2,ind) = t;
    end
end
shaperConv = sortrows(shaperConv',2)'
