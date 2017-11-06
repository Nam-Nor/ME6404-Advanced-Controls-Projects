% crane trajectory creation
% Written by my lab partner, Lucas Tiziani
clear, clc, clf

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The trajectory waypoints are mapped out below:

% Trajectory
% 1. [65,567,1720] for start location
% 2. 821 (mm) max hoist position for obstacle clearance, set to 811
% 3. [307,868,811] for payload r/theta location
% 4. 1561 (mm) min hoist pos for magnet engage, set to 1565
% 5. 879 (mm) max hoist poisition for payload pickup off ground, set to 835
% 6. max trolley position for obstacle avoidance, TODO: check max
% 7. [197,-,835] for payload dropoff location
% 8. 811 (mm) trolley location for payload dropoff
% 9. 988 mm hoist position for payload dropoff

% trajectory waypoints
%          1     2      3   4      5    6    7    8    9
posDes = [[65,   -1,  307,  -1,    -1,  -1,  197, -1,  -1];        % slew (deg)
          [567,  -1,  868,  -1,    -1,  900, -1,  811, -1]/1000;   % rad (m)
          [1720, 811, -1,   1580,  835, -1,  835, -1,  988]/1000]; % hoist (m)
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      
% crane parameters
tSample = 0.04; % (s) sample time

% max velocities
velMax = [0.35*(180/pi); % slew (deg/s)
          0.14;          % rad (m/s) 
          0.2];         % hoist (m/s) TODO: determine this value
      
% max accelerations
accMax = [0.7*(180/pi); % slew (deg/s^2)
          1.20;         % rad (m/s^2)
          1.2];        % hoist (m/s^2) TODO: determine this value

% calculate acceleration/deceleration time & distance
tVelRamp = velMax./accMax; % (s)
posVelRamp = (1/2)*accMax.*tVelRamp.^2; % (deg)

% initialize trajectory arrays
tTraj = zeros(length(posDes),4,3);
vTraj = zeros(length(posDes),4,3);
aTraj = zeros(length(posDes),4,3);
dist = zeros(3,1); % travel distances
dir = zeros(3,1); % travel directions


%% calculate accelerations and corresponding times based on waypoints
pos = posDes(:,1); % set initial positions to first waypoint
for i = 2:length(posDes)
%     fprintf('%i\n',i)
    tMax = max(tTraj(i-1,end,:)); % end time of last waypoint
    
% calculate travel distance and direction
    for j = 1:3
        if posDes(j,i) < 0 % if no change to desired location
            dist(j) = 0; % set distance to zero
        else % otherwise
            dist(j) = posDes(j,i) - pos(j); % (deg) set distance
            dir(j) = sign(dist(j)); % determine direction based on sign
        end
    
        if dist(j) == 0 % if no distance to travel
            tTraj(i,:,j) = tMax; % just repeat time from last timestep
        elseif abs(dist(j)) < 2*posVelRamp(j) % if max velocity not reached
            t = sqrt(abs(dist(j))/accMax(j));
            tTraj(i,:,j) = [0,t,t,2*t] + tMax;
            v = dir(j)*accMax(j)*t;
            vTraj(i,:,j) = [0,v,v,0];
            aTraj(i,:,j) = [dir(j)*accMax(j), 0, -dir(j)*accMax(j), 0];
        else % if max velocity reached
            t = (abs(dist(j)) - 2*posVelRamp(j))/velMax(j);
            tTraj(i,:,j) = [0, tVelRamp(j), tVelRamp(j) + t,... 
                t + 2*tVelRamp(j)] + tMax;
            vTraj(i,:,j) = [0, dir(j)*velMax(j), dir(j)*velMax(j), 0];
            aTraj(i,:,j) = [dir(j)*accMax(j), 0, -dir(j)*accMax(j), 0];
           
        end
        
        % update positions to current waypoint
        if posDes(j,i) >= 0
            pos(j) = posDes(j,i);
        end
        
    end
end

% reshape data into vectors
tSlew = reshape(tTraj(:,:,1)',[],1);
vSlew = reshape(vTraj(:,:,1)',[],1);
aSlew = reshape(aTraj(:,:,1)',[],1); 
tRad = reshape(tTraj(:,:,2)',[],1);
aRad = reshape(aTraj(:,:,2)',[],1);
tHoist = reshape(tTraj(:,:,3)',[],1);
aHoist = reshape(aTraj(:,:,3)',[],1);

% group vectors into cell array
tArray = {tSlew,tRad,tHoist};
aArray = {aSlew,aRad,aHoist};

% eliminate repeated times
for j = 1:3
    tVecRep = tArray{j}; % vector of times with repeats
    aVecRep = aArray{j}; % vector of corresponding accel values
    
    [~,~,Iu] = unique(tVecRep); % get indices corresponding to repeated times
    tVec = zeros(Iu(end),1); % initialize vector of non-repeated times
    aVec = zeros(Iu(end),1); % initialize corresponding accel vector


    for i = 1:Iu(end) % loop over all unique time values
        
        inds = (Iu==i); % indices of repeated time values
        c = sum(inds); % number of repeated time values
        if (c > 1) % if multiple occurrences
            tRep = tVecRep(inds); % repeated time values
            aRep = aVecRep(inds); % corresponding accel values

            [~,I] = max(abs(aVecRep(Iu==i))); % index of accel value to keep

            tVec(i,1) = tRep(I);
            aVec(i,1) = aRep(I);
        else
            tVec(i,1) = tVecRep(inds);
            aVec(i,1) = aVecRep(inds);
        end
    end
    tArr{j} = tVec;
    aArr{j} = aVec;
end


%% calculate acceleration, velocity, position arrays
% loop over sample times
tInc = 0.0001; % small time increment
t = 0:tInc:(max(tTraj(end,4,:)) + 0.1);
a = zeros(length(t),3);
v = zeros(length(t),3);
p = zeros(length(t),3);
p(1,:) = posDes(:,1)';
ind = [1,1,1]; % index of current accel

for i = 2:(length(t))
    tCur = t(i);
    
    for j = 1:3 % loop over 3 crane axes
        if ind(j) < length(tArr{j}) % if final time not yet reached
            if tCur >= tArr{j}(ind(j)+1) % if next time reached
                ind(j) = ind(j) + 1; % increment index          
            end
        end
        
        a(i,j) = aArr{j}(ind(j));
        v(i,j) = a(i,j)*tInc + v(i-1,j);
        p(i,j) = v(i,j)*tInc + p(i-1,j);
    end   
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% at this point in the script, arrays a,v, and p should have fairly
% accurate data for unshaped trajectories (the timestep for these
% trajectories is 'tInc'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% shape velocity TODO: adjust for max velocity not reached
w1 = sqrt(9.81/0.821); % natural frequency
T1 = 2*pi/w1; % period
w2 = sqrt(9.81/0.7); % natural frequency
T2 = 2*pi/w2; % period

tShaper = 0:tInc:(T1/2); % shaper time vector
ZV = [0.5; zeros(length(tShaper)-2,1); 0.5]; % shaper values
velShaped{1} = conv(v(:,1),ZV); % shaped velocity
velShaped{2} = conv(v(:,2),ZV); 

timeShaped{1} = 0:tInc:(length(velShaped{1})*tInc - tInc);
timeShaped{2} = 0:tInc:(length(velShaped{2})*tInc - tInc);


% velShaped = {};
% timeShaped = {};
% for j = 1:2 % loop over slew and rad
%     % determine nonzero velocity command times
%     i = 1;
%     velTimes = [];
%     while i <= length(tArr{j})
%         if aArr{j}(i) ~= 0
%             velTimes = [velTimes, tArr{j}(i), tArr{j}(i+3)];
%             i = i + 4;
%         else
%             i = i + 1;
%         end
%     end
% 
% % shape each velocity command separately based on hoist length
%     vShaped = v(:,j);
%     tShaped = t;
%     tAdded = 0;
%     for i = 1:2:length(velTimes)
%         tStart = velTimes(i) + tAdded; % velocity command start time
%         tEnd = velTimes(i+1) + tAdded; % velocity command end time
% 
%         [~,indStart] = min(abs(t-tStart)); % corresponding start index
%         [~,indEnd] = min(abs(t-tEnd)); % corresponding end index
% 
%         w = sqrt(9.81/p(indStart,3)); % natural frequency
%         T = 2*pi/w; % period
% 
%         tShaper = 0:tInc:(T/2); % shaper time vector
%         ZV = [0.5; zeros(length(tShaper)-2,1); 0.5]; % shaper values
%         vS = conv(vShaped(indStart:indEnd,1),ZV); % shaped velocity
% 
%         vShaped = [vShaped(1:(indStart-1)); vS(:,1); vShaped((indEnd+1):end)]; % new velocity vector
%         tShaped = 0:tInc:(length(vShaped)*tInc - tInc); % new time vector
%         
%         tAdded = tShaped(end) - t(end); % time added by shaper convolution
%     end
%     velShaped{j} = vShaped;
%     timeShaped{j} = tShaped;
% end


%% discretize to crane sample rate
t = t(1:(tSample/tInc):end);
a = a(1:(tSample/tInc):end,:);
v = v(1:(tSample/tInc):end,:);
p = p(1:(tSample/tInc):end,:);


%% plot
% unshaped vs shaped velocity
subplot(3,1,1)
plot(t,v(:,1))
hold on
plot(timeShaped{1},velShaped{1})
grid on
legend('unshaped','shaped')
title('Slew')

subplot(3,1,2)
plot(t,v(:,2))
hold on
plot(timeShaped{2},velShaped{2})
grid on
legend('unshaped','shaped')

subplot(3,1,3)
plot(t,p(:,3))
hold on
plot(t,v(:,3))
plot(t,a(:,3))
grid on
xlim([0 max(t)])
legend('position','velocity','accleration')
title('Hoist')


% position, velocity, acceleration
% subplot(3,2,1)
% plot(t,p(:,1))
% hold on
% plot(t,v(:,1))
% plot(t,a(:,1))
% grid on
% xlim([0 max(t)])
% legend('position','velocity','accleration')
% title('Slew')
% 
% subplot(3,2,3)
% plot(t,p(:,2))
% hold on
% plot(t,v(:,2))
% plot(t,a(:,2))
% grid on
% xlim([0 max(t)])
% legend('position','velocity','accleration')
% title('Radius')
% 
% subplot(3,2,5)
% plot(t,p(:,3))
% hold on
% plot(t,v(:,3))
% plot(t,a(:,3))
% grid on
% xlim([0 max(t)])
% legend('position','velocity','accleration')
% title('Hoist')
% 
% subplot(3,2,[2,4])
% polarplot(p(:,1)*pi/180,p(:,2))


v(:,1) = v(:,1)/velMax(1);
v(:,2) = v(:,2)/velMax(2);
v(:,3) = v(:,3)/velMax(3)*-1;
v = round(v*100);

