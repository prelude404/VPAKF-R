function [x_esti] = kf_sa(gtd, t, vel, uwb)
% kf_sa: The classical Kalman filter with state augmentation
% System Model:
%   x(k) = A x(k-1) + B u(k) + q
%   y(k) = H x(k) + r
% State Augmentation:
%   z = [P', P0'B, ||B||^2, B']

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
A = [eye(3), zeros(3,2), -dt*eye(3); zeros(5,3), eye(5)];

B = [dt*eye(3); zeros(5,3)];

%% Initialization
bias = zeros(3,K);
% % bag1
% bias(:,1) = [0.053;-0.012;0.020];

% bag2
bias(:,1) = [-0.043;0.064;-0.124];

% % bag5
% bias(:,1) = [0.031;-0.067;0.100];

% bias(:,1) = [0;0;0];
x0 = [gtd(1:3,1)', gtd(1:3,1)'*bias(:,1), norm(bias(:,1))^2, bias(:,1)']';
x_esti = zeros(8,K);
x_predict = zeros(8,K);
x_pre = zeros(8,K);
x_esti(:,1) = x0;
% x_predict(:,1) = x_esti(:,1);
x_pre(:,1) = x_esti(:,1);

%% Parameter
% % bag1
% R = 1e-2;
% P = diag([1e-5*[1,1,1],1e-1,1e-8,1e-7*[1,1,1]]);
% Q = B * diag(1e-3*[0.01,0.01,0.01]) * B';

% bag2
R = 1e-2;
P = diag([1e-3*[1,1,1],1e0,1e-9,1e-7*[1,1,1]]);
Q = B * diag(1e-7*[0.01,0.01,0.01]) * B';

% % bag5
% R = 1e-2;
% P = diag([1e-5*[1,1,1],1e0,1e-10,1e-7*[1,1,1]]);
% Q = B * diag(1e-4*[0.01,0.01,0.01]) * B';

lambda = 1;

%% KF_SA
for i = 2:K
    x_pre(:,i) = A * x_esti(:,i-1) + B * v(:,i);
    P = A * P * A'+Q;
    H = [2*Iv(1,i), 2*Iv(2,i), 2*Iv(3,i), -2*(i-1)*dt, ((i-1)*dt)^2, 0,0,0];
    K = P * H' * ((H*P*H'+R)^(-1));
    x_esti(:,i) = x_pre(:,i) + K * (y(i) - H * x_pre(:,i));
    P = 1/lambda * (eye(8)-K*H)*P;
end


