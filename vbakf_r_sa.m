function [x_esti] = vbakf_r_sa(gtd, t, vel, uwb)
% vbakf_q: The variational Bayesian adaptive Kalman filter with unknown R.
% R: covariance of measurement noise
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
bias(:,1) = [0.0374;-0.0088;0.0063];
% bias(:,1) = [0;0;0];
x0 = [gtd(1:3,1)', gtd(1:3,1)'*bias(:,1), norm(bias(:,1))^2, bias(:,1)']';
x_esti = zeros(8,K);
x_pre = zeros(8,K);
x_esti(:,1) = x0;
x_pre(:,1) = x_esti(:,1);

%% Parameter
alpha =1e4;
beta = 1e2;
rou = 0.9999;
N = 10;
P = diag([1e-5*[1,1,1],1e-1,1e-8,1e-7*[1,1,1]]);
Q = B * diag(1e-5*[1,1,1]) * B';

%% VBAKF-R
for i = 2:K
    H = [2*Iv(1,i), 2*Iv(2,i), 2*Iv(3,i), -2*(i-1)*dt, ((i-1)*dt)^2, 0,0,0];
    
    % Prediction
    x_pre(:,i) = A * x_esti(:,i-1) + B * v(:,i);
    P_pre = A * P * A' + Q;
    
    % Initialization
    P = P_pre;
    alpha = rou * alpha + 0.5;
    beta_pre = rou * beta;
    beta = beta_pre;
    for j = 1:N-1
        sigma = beta / alpha;
        x = x_pre(:,i) + P_pre * H' / (H * P_pre * H' + sigma) * (y(i) - H * x_pre(:,i));
        P = P_pre - P_pre * H' / (H * P * H' + sigma) * H * P_pre;
        beta = beta_pre + 0.5 * (y(i) - H * x)^2 + 0.5 * (H * P * H');
    end
    R(i) = beta / alpha;
    x_esti(:,i) = x;
end

figure(11)
plot(t,R,'b-','linewidth',1);
legend('R','FontName','Times New Roman','FontSize',12);
set(gca,'linewidth',0.5,'FontSize',8,'FontName','Times New Roman');
ylabel('Measurement Noise Covariance (m^2)','FontName','Times New Roman','FontSize',16);
xlabel('Time (s)','FontName','Times New Roman','FontSize',16);
grid on
% title('Measurement Noise Covariance','FontName','Times New Roman','FontSize',16);





