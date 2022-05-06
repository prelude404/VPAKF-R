function [x_esti,x_predict] = kf(gtd, t, vel, uwb)
% kf: The classical Kalman filter.
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
R = zeros(1,K);

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
x_predict = zeros(3,K);
x_pre = zeros(3,K);
x_esti(:,1) = x0;
x_predict(:,1) = x_esti(:,1);
x_pre(:,1) = x_esti(:,1);

%% Parameter
% bag1
R = 1e-1;
P = diag(1e-5*[1,1,1]);
Q = B * diag(1e-5*[1,1,1]) * B';

% % bag2
% R = 1e-4;
% P = diag(1e-4*[1,1,1]);
% Q = B * diag(1e-5*[1,1,1]) * B';

%% Bias
bias = zeros(3,K);
bias(:,1) = [0.0374;-0.0088;0.0063];
x_bias = zeros(3,K);
% should find a physical explanation for a and b
a = 1-1e-2;
b = 8;

%% KF
for i = 2:K
%     v(:,i) = v(:,i) - bias(:,i-1);
    x_predict(:,i) = A * x_predict(:,i-1) + B * v(:,i);
    x_pre(:,i) = A * x_esti(:,i-1) + B * v(:,i);
    P = A * P * A'+Q;
    H = [Iv(1,i),Iv(2,i),Iv(3,i)];
    K = P * H' * ((H*P*H'+R)^(-1));
    x_esti(:,i) = x_pre(:,i) + K * (y(i) - H * x_pre(:,i));
    P = (eye(3)-K*H)*P;
    x_bias(:,i) = x_pre(:,i) - x_esti(:,i);
    bias(:,i) = a * bias(:,i-1) + b * x_bias(:,i);
end
% figure(10)
% plot(t,bias(1,:),'-r',t,bias(2,:),'g-',t,bias(3,:),'b-');