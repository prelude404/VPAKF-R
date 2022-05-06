function [x_esti] = vbakf_pm(gtd, t, vel, uwb)
% vbakf_pm: The variational Bayesian adaptive Kalman filter with unknown
% process and measurement noise covariance matrices.
% P: process noise covariance
% M: measurement noise covariance
% System Model:
%   x(k) = A x(k-1) + B u(k) + q
%   y(k) = H x(k) + r

%% Data Preparation: y,v,K,dt
x0 = gtd(1:3,1);
K = length(t);
dt = t(2) - t(1);
v = vel;
Iv = zeros(3,K);
y = zeros(1,K);
y0 = uwb(1).^2;
y1 = uwb.*uwb;

for i=1:K
    Iv(1,i) = dt * trapz(v(1,1:i));
    Iv(2,i) = dt * trapz(v(2,1:i));
    Iv(3,i) = dt * trapz(v(3,1:i));
    
    y(i) = 0.5 * y1(i) - 0.5 * y0 + 0.5 * (Iv(1,i)^2 + Iv(2,i)^2 + Iv(3,i)^2);
end

%% x_k = A x_k-1 + B v_k
A = [1,0,0;
     0,1,0;
     0,0,1];

B = [dt,0,0;
     0,dt,0;
     0,0,dt];

%% Initialization
x0 = x0 + [0;0;0]; % add initial error
x_esti = zeros(3,K);
x_pre = zeros(3,K);
x_esti(:,1) = x0;
x_pre(:,1) = x_esti(:,1);

%% Parameter
P = 1e-5 * eye(3);
Q = B * (1e-3*eye(3)) * B';
R0 = 1e-2;
tau = 3;
rou = 1-1e-2;
u = 2.1;

n = 3;
m = 1;
U = (u-m-1) * R0;
N = 10;

%% VBAKF-PM
for i = 2:K
%     v(:,i) = v(:,i)-[0.001;-0.002;0.0005];
%     v(:,i) = v(:,i) - bias(:,i-1);
    % Prediction
    x_pre(:,i) = A * x_esti(:,i-1) + B * v(:,i);
    P_pre = A * P * A' + Q;
    
    x = x_pre(:,i);
    P = P_pre;
    t_pre = n + tau + 1;
    T_pre = tau * P_pre;
    u_pre = rou * (u-m-1) + (m+1);
    U_pre = rou * U;
    
    H = [Iv(1,i),Iv(2,i),Iv(3,i)];
    
    for j = 1:N
        % update Pk|k-1
        A = P + (x-x_pre(:,i)) * (x-x_pre(:,i))';
        t = t_pre + 1;
        T = A + T_pre;
        % update Rk
        B = (y(i)-H*x)*(y(i)-H*x)' + H*P*H';
        u = u_pre +1;
        U = B + U_pre;
        % update xk
        P = ((t-n-1)*(T^(-1)))^(-1);
        R = ((u-m-1)*(U^(-1)))^(-1);
        K = P * H' * ((H*P*H'+R)^(-1));
        x = x_pre(:,i) + K * (y(i)-H*x_pre(:,i));
        P = P - K*H*P;
    end
    
    x_esti(:,i) = x;
    
%     x_bias(:,i) = x_pre(:,i) - x_esti(:,i);
%     bias(:,i) = a * bias(:,i-1) + b * x_bias(:,i);
end





