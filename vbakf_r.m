function [x_esti,x_predict] = vbakf_r(gtd, t, vel, uwb)
% vbakf_q: The variational Bayesian adaptive Kalman filter with unknown R.
% R: covariance of measurement noise
% System Model:
%   x(k) = A x(k-1) + B u(k) + q
%   y(k) = H x(k) + r

%% Data Preparation: v,K,dt
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
alpha =1e4;
beta = 1e3;
rou = 0.99;
N = 10;
P = diag(1e-5*[1,1,1]);
% Q = diag([1e-7*[1,1,1],1e-6*[1,1,1],1e-3,1e-7]);
Q = B * diag(1e-3*[0.01,0.01,0.01]) * B';

%% Bias
bias = zeros(3,K);
bias(:,1) = [0.001;-0.002;0.0005];
x_bias = zeros(3,K);
% should find a physical explanation for a and b
a = 0.92;
b = 3;
% a = 0.92;
% b = 10;

%% VBAKF-R
for i = 2:K
%     v(:,i) = v(:,i)-[0.001;-0.002;0.0005];
    x_predict(:,i) = A * x_predict(:,i-1) + B * v(:,i);
    % Prediction
    v(:,i) = v(:,i) - bias(:,i-1);
    x_pre(:,i) = A * x_esti(:,i-1) + B * v(:,i);
    P_pre = A * P * A' + Q;
    P = P_pre;
    alpha = rou * alpha + 0.5;
    beta_pre = rou * beta;
    beta = beta_pre;
    H = [Iv(1,i),Iv(2,i),Iv(3,i)];
    for j = 1:N-1
        sigma = beta / alpha;
        x = x_pre(:,i) + P_pre * H' / (H * P_pre * H' + sigma) * (y(i) - H * x_pre(:,i));
        P = P_pre - P_pre * H' / (H * P * H' + sigma) * H * P_pre;
        beta = beta_pre + 0.5 * (y(i) - H * x)^2 + 0.5 * (H * P * H');
    end
    R(i) = beta / alpha;
    x_esti(:,i) = x;
    x_bias(:,i) = x_pre(:,i) - x_esti(:,i);
    bias(:,i) = a * bias(:,i-1) + b * x_bias(:,i);
end

% figure(10)
% plot(t,bias,'linewidth',1);
% legend('bias_X','bias_Y','bias_Z','FontName','Times New Roman','FontSize',12);
% ylabel('Bias','FontName','Times New Roman','FontSize',12);
% xlabel('Time','FontName','Times New Roman','FontSize',12);
% title('Bias Estimation','FontName','Times New Roman','FontSize',16);
% disp(['X: ',num2str(mean(bias(1,:))),' Y: ',num2str(mean(bias(2,:))),' Z: ',num2str(mean(bias(3,:)))]);

% figure(11)
% plot(t,R,'linewidth',1);
% legend('cov_R','FontName','Times New Roman','FontSize',12);
% ylabel('Covariance','FontName','Times New Roman','FontSize',12);
% xlabel('Time','FontName','Times New Roman','FontSize',12);
% title('Measurement Noise Covariance','FontName','Times New Roman','FontSize',16);






