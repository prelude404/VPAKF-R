function [x_esti] = r_sa_sw(gtd, t, vel, uwb)
% vbakf_r: The variational Bayesian adaptive Kalman filter with unknown R.
% R: covariance of measurement noise
% System Model:
%   x(k) = A x(k-1) + B u(k) + q
%   y(k) = H x(k) + r
% State Augmentation:
%   z = [P', P0'B, ||B||^2, B']
% Sliding Window:
% periodically reset the initial measurement y0

%% Data Preparation: v,K,dt
K = length(t);
dt = t(2) - t(1);
v = vel;
Iv = zeros(3,K);
y = zeros(1,K);
R = zeros(1,K);

win_len = 25;

%% Initialization
bias = zeros(3,K);
x_esti = zeros(8,K);
x_pre = zeros(8,K);
for i=1:win_len
    bias(:,i) = [0.053;0.000;0.010];
    x_esti(:,i) = [gtd(1:3,i)', gtd(1:3,1)'*bias(:,i), norm(bias(:,i))^2, bias(:,i)']';
end
x_pre(:,1:win_len) = x_esti(:,1:win_len);

%% x_k = A x_k-1 + B v_k
A = [eye(3), zeros(3,2), -dt*eye(3); zeros(5,3), eye(5)];
B = [dt*eye(3); zeros(5,3)];

%% Parameter
alpha = 1e4;
beta = 1e2;
rou = 0.99;
N = 10;

P = diag([1e-4*[1,1,1],1e-7,1e-7,1e-8*[1,1,1]]);
Q = B * diag(1e-3*[0.01,0.01,0.01]) * B';

%% KF_SA_SW
for i = (win_len+1):K
    if rem(i,win_len)==1
        % reset
        start_i = i-1;
        p0 = x_esti(1:3,i-1);
        y0 = uwb(i-1).^2;
        x_esti(4,i-1) = p0' * x_esti(6:8,i-1);
    end
    Iv(1,i) = dt * trapz(v(1,start_i:i));
    Iv(2,i) = dt * trapz(v(2,start_i:i));
    Iv(3,i) = dt * trapz(v(3,start_i:i));
    y(i) = uwb(i).*uwb(i) - y0 + (Iv(1,i)^2 + Iv(2,i)^2 + Iv(3,i)^2);
    % predict
    x_pre(:,i) = A * x_esti(:,i-1) + B * v(:,i);
    P_pre = A * P * A'+Q;
    % update
    P = P_pre;
    alpha = rou * alpha + 0.5;
    beta_pre = rou * beta;
    beta = beta_pre;
    H = [2*Iv(1,i), 2*Iv(2,i), 2*Iv(3,i), -2*(i-start_i)*dt, ((i-start_i)*dt)^2, 0,0,0];
    for j = N-1
        sigma = beta / alpha;
        x = x_pre(:,i) + P_pre * H' / (H * P_pre * H' + sigma) * (y(i) - H * x_pre(:,i));
        P = P_pre - P_pre * H' / (H * P * H' + sigma) * H * P_pre;
        beta = beta_pre + 0.5 * (y(i) - H * x)^2 + 0.5 * (H * P * H');
    end
    R(i) = beta / alpha;
    x_esti(:,i) = x;
end

