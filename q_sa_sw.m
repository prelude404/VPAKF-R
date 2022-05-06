function [x_esti] = q_sa_sw(gtd, t, vel, uwb)
% vbakf_q: The variational Bayesian adaptive Kalman filter with unknown Q.
% Q: covariance of process noise
% System Model:
%   x(k) = A x(k-1) + B u(k) + q
%   y(k) = H x(k) + r
% State Augmentation:
%   z = [P', P0'B, ||B||^2, B']
% Sliding Window:
% periodically reset the initial measurement y0

%% Data Preparation: y,v,K,dt
K = length(t);
dt = t(2) - t(1);
v = vel;
Iv = zeros(3,K);
y = zeros(1,K);

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
rho = 0.99;
nx = 8;
tau = 20;

R = 1e-2;
P = diag([1e-6*[1,1,1],1e-7,1e-6,1e-8*[1,1,1]]);
Q0 = diag([1e-7*[1,1,1],1e-7,1e-7,1e-8*[1,1,1]]);

mu = tau + nx +1;
U = tau * Q0;
delta = 1e-6;
N = 3;


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
    sigma = A * P * A';
    mu = rho * (mu-nx-1) + nx + 1;
    U = rho * U;
    % update
    H = [2*Iv(1,i), 2*Iv(2,i), 2*Iv(3,i), -2*(i-start_i)*dt, ((i-start_i)*dt)^2, 0,0,0];
    theta = x_pre(:,i);
    x_former = theta;
    
    for j = 1:N
        Aq = U./mu;
        % update x
        K = Aq * H' * ((H*Aq*H'+R)^(-1));
        x = theta + K * (y(:,i) - H * theta);
        P = Aq - K * H * Aq;
        % update theta
        K1 = sigma * ((sigma + Aq)^(-1));
        theta = x_pre(:,i) + K1 * (x - x_pre(:,i));
        P1 = sigma - K1 * sigma;
        % update Q
        mu = mu + 1;
        U = U + (x-theta)*(x-theta)' + P + P1;
        if norm(x - x_former)/norm(x) < delta
            break;
        end
        x_former = x;
    end
    x_esti(:,i) = x;
end


