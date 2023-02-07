clear all
close all

%importing torque values
u = importdata('torques.txt');

%setting initial conditions and q arrays
IC1 = [0;0;0;0];
IC2 = [0;0;0;0];
q1 = zeros(length(u),2);
q2 = zeros(length(u),2);

%options variable for ode45
options = odeset('RelTol',1e-4,'AbsTol',[1e-5 1e-5 1e-5 1e-5]);

for i = 1:length(u)-1
    %1 kHz
    [tsim,y1] = ode45(@samplesys,[0 1/1000],IC1,options,u(i,:)');
    IC1 = y1(end,:);
    q1(i+1,:) = y1(end,1:2);
    
    %200 Hz
    [tsim,y2] = ode45(@samplesys,[0 1/200],IC2,options,u(i,:)');
    IC2 = y2(end,:);
    q2(i+1,:) = y2(end,1:2);
end

%plotting both figures
figure(1)
pos1 = fwd_kin(q1(:,1),q1(:,2));
plot(pos1(:,1),pos1(:,2));
title('Part 1: 1 kHz');
xlabel('X [m]');
ylabel('Y [m]');

figure(2)
pos2 = fwd_kin(q2(:,1),q2(:,2));
plot(pos2(:,1),pos2(:,2));
title('Part 2: 200 Hz');
xlabel('X [m]');
ylabel('Y [m]');

%state space function
function dx = samplesys(t,x,T)
%taking x = [q1; q2; q_1dot; q_2dot]

%calling dynamics_matrics for BCG
[B, C, G] = dynamics_matrices([x(1);x(2)],[x(3);x(4)]);

%getting qdoubledot to get dx(2) and dx(4)
qdoubledot = B\(T-G-C*[x(3);x(4)]);
dx = [x(3);x(4); qdoubledot];
end
