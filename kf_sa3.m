function [x_esti] = kf_sa3(gtd, t, vel, uwb)
% kf_sa: The classical Kalman filter with state augmentation
% System Model:
%   x(k) = A x(k-1) + B u(k) + q
%   y(k) = H x(k) + r
% State Augmentation:
%   x = [r', vf', y, r'vf, ||vf||^2]'

%% Data Preparation: y,v,K,dt
K = length(t);
dt = t(2) - t(1);
v = vel;
y = uwb;

%% x_k = A x_k-1 + B v_k
B = [diag([dt,dt,dt]);zeros(6,3)];
%% y_k = H x_k
H = [0,0,0,0,0,0,1,0,0];

%% Initialization
bias = zeros(3,K);

% % bag1
% bias(:,1) = [0.043;-0.012;0.031];

% bag2
bias(:,1) = [-0.204;0.000;-0.196];

% % bag5
% bias(:,1) = [0.031;-0.067;0.100];

x0 = [gtd(1:3,1)', bias(:,1)', norm(gtd(1:3,1)), gtd(1:3,1)'*bias(:,1), bias(:,1)'*bias(:,1)]';
x_esti = zeros(9,K);
x_pre = zeros(9,K);
x_esti(:,1) = x0;
x_pre(:,1) = x0;

%% Parameter
% % bag1
% R = 1e1;
% Q = B * diag(1e-6*[1,1,1]) * B';
% P = diag([1e-7*[1,1,1],1e-7*[1,1,1],1e-5,1e-2,1e-9]);

% bag2
R = 1e-3;
Q = B * diag(1e-3*[1,1,1]) * B';
P = diag([1e-3*[1,1,1],1e-5*[1,1,1],1e0,1e0,1e-11]);

% bag5

lambda = 1;
%% KF_SA
for i = 2:K
    A = [eye(3),-dt*eye(3),zeros(3,3);
         zeros(3,3),eye(3),zeros(3,3);
         dt/y(i)*v(:,i)',zeros(1,3),1,-dt/y(i),0;
         zeros(1,3),v(:,i)',0,1,-dt;
         zeros(1,8),1];
    x_pre(:,i) = A * x_esti(:,i-1) + B * v(:,i);
    P = A * P * A'+Q;
    K = P * H' * ((H*P*H'+R)^(-1));
    x_esti(:,i) = x_pre(:,i) + K * (y(i) - H * x_pre(:,i));
    P = 1/lambda * (eye(9)-K*H)*P;
end


