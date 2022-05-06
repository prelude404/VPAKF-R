function [x_esti] = kf_sa_sw(gtd, t, vel, uwb)
% kf_sa: The classical Kalman filter with state augmentation
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
%     bias(:,i) = [0.053;0.000;0.010];
    bias(:,i) = [0;0;0];
    x_esti(:,i) = [gtd(1:3,i)', gtd(1:3,1)'*bias(:,i), norm(bias(:,i))^2, bias(:,i)']';
end
x_pre(:,1:win_len) = x_esti(:,1:win_len);

%% x_k = A x_k-1 + B v_k
A = [eye(3), zeros(3,2), -dt*eye(3); zeros(5,3), eye(5)];
B = [dt*eye(3); zeros(5,3)];

%% Parameter
R = 0.05;
P = diag([1e-4*[1,1,1],1e-7,1e-7,1e-8*[1,1,1]]);
Q = B * diag(1e-3*[0.01,0.01,0.01]) * B';
lambda = 1;

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
    
    x_pre(:,i) = A * x_esti(:,i-1) + B * v(:,i);
    P = A * P * A'+Q;
    H = [2*Iv(1,i), 2*Iv(2,i), 2*Iv(3,i), -2*(i-start_i)*dt, ((i-start_i)*dt)^2, 0,0,0];
    K = 1/lambda * P * H' * ((H*P*H'+R)^(-1));
    x_esti(:,i) = x_pre(:,i) + K * (y(i) - H * x_pre(:,i));
    P = (eye(8)-K*H)*P;
end

% figure(10)
% plot(t,x_esti(6,:),'-r',t,x_esti(7,:),'g-',t,x_esti(8,:),'b-');

        






