function [x_esti] = vbakf_q(gtd, t, vel, uwb)
% vbakf_q: The variational Bayesian adaptive Kalman filter with unknown Q.
% Q: covariance of process noise
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
rho = 0.99;
nx = 3;
tau = 0.1;

R = 1e-2;
P = diag(1e-5*[1,1,1]);
Q0 = B * diag(1e-3*[0.01,0.01,0.01]) * B';

mu = tau + nx +1;
U = tau * Q0;
delta = 1e-6;
N = 5;

%% Bias
bias = zeros(3,K);
bias(:,1) = [0.0374;-0.0088;0.0063];
x_bias = zeros(3,K);
% should find a physical explanation for a and b
a = 0.93;
b = 3;

%% VBAKF-Q
for i = 2:K
    % Prediction
    v(:,i) = v(:,i) - bias(:,i-1);
    x_pre(:,i) = A * x_esti(:,i-1) + B * v(:,i);
    sigma = A * P * A';
    mu = rho * (mu-nx-1) + nx + 1;
    U = rho * U;
    % Update
    H = [Iv(1,i),Iv(2,i),Iv(3,i)];
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
    x_bias(:,i) = x_pre(:,i) - x_esti(:,i);
    bias(:,i) = a * bias(:,i-1) + b * x_bias(:,i);
end




