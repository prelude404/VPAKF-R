function [x_esti] = kf_sa2(gtd, t, vel, uwb)
% kf_sa: The classical Kalman filter with state augmentation
% System Model:
%   x(k) = A x(k-1) + B u(k) + q
%   y(k) = H x(k) + r
% State Augmentation:
%   z = [P', P0'b0, P0'k, b0'k, ||b0||^2, ||k||^2, b0', k']

%% Data Preparation: y,v,K,dt   
K = length(t);
dt = t(2) - t(1);
v = vel;
Iv = zeros(3,K);
y = zeros(1,K);
y0 = uwb(1).^2;
y1 = uwb.*uwb;
R = zeros(1,K);

for i=1:K
    Iv(1,i) = dt * trapz(v(1,1:i));
    Iv(2,i) = dt * trapz(v(2,1:i));
    Iv(3,i) = dt * trapz(v(3,1:i));
    
    y(i) = y1(i) - y0 + (Iv(1,i)^2 + Iv(2,i)^2 + Iv(3,i)^2);
end

%% x_k = A x_k-1 + B v_k
% A = [eye(3), zeros(3,2), -dt*eye(3); zeros(5,3), eye(5)];

B = [dt*eye(3); zeros(11,3)];

%% Initialization
bias = zeros(3,K);
b0 = zeros(3,K);
k = zeros(3,K);

% bag1
bias(:,1) = [0.050;-0.012;0.021];
k(:,1) = [3e-4;-1e-4;-1e-4];

% % bag2
% bias(:,1) = [-0.016;0.069;-0.123];
% k(:,1) = [3e-4;1e-4;-12e-4];

b0(:,1) = bias(:,1);

x0 = [gtd(1:3,1)', gtd(1:3,1)'*b0(:,1), gtd(1:3,1)'*k(:,1), b0(:,1)'*k(:,1), norm(b0(:,1))^2, norm(k(:,1))^2, b0(:,1)', k(:,1)']';
x_esti = zeros(14,K);
x_pre = zeros(14,K);
x_esti(:,1) = x0;
x_pre(:,1) = x_esti(:,1);

%% Parameter
% bag1
R = 1e1;
P = diag([1e-2*[1,1,1],1e1,1e-5,1e-6,1e-4,1e-3,1e-4*[1,1,1],1e-11*[1,1,1]]);
Q = B * diag(1e-5*[0.01,0.01,0.01]) * B';

% % bag2
% R = 1e1;
% P = diag([1e-1*[1,1,1],1e2,1e-5,1e-7,1e-7,1e-3,1e-4*[1,1,1],1e-8*[1,1,1]]);
% Q = B * diag(1e-0*[0.01,0.01,0.01]) * B';

lambda = 1;

%% KF_SA
for i = 2:K
    A = [eye(3), zeros(3,5), -dt*eye(3), -dt*(i-1)*dt*eye(3); zeros(11,3), eye(11)];
    H = [2*Iv(:,i)',-2*(i-1)*dt,-((i-1)*dt)^2,((i-1)*dt)^3,((i-1)*dt)^2,0.25*((i-1)*dt)^2,zeros(1,6)];
    
    x_pre(:,i) = A * x_esti(:,i-1) + B * v(:,i);
    P = A * P * A'+Q;
    K = P * H' * ((H*P*H'+R)^(-1));
    x_esti(:,i) = x_pre(:,i) + K * (y(i) - H * x_pre(:,i));
    P = 1/lambda * (eye(14)-K*H)*P;
    
    bias(:,i) = x_esti(9:11,i) + (i-1)*dt*x_esti(12:14,i);
end

figure(11)
plot(t,bias(1,:),'-r',t,bias(2,:),'g-',t,bias(3,:),'b-','linewidth',1);
title('Estimation Bias','FontName','Times New Roman','FontSize',16);
xlabel('Time (s)','FontName','Times New Roman','FontSize',12);
ylabel('Error (m/s)','FontName','Times New Roman','FontSize',12);
legend('X_{bias}','Y_{bias}','Z_{bias}','FontName','Times New Roman','FontSize',8);
