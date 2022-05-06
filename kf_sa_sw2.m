function [x_esti] = kf_sa_sw2(gtd, t, vel, uwb)
% kf_sa: The classical Kalman filter with state augmentation
% System Model:
%   x(k) = A x(k-1) + B u(k) + q
%   y(k) = H x(k) + r
% State Augmentation:
%   z = [P', P0'b0, P0'k, b0'k, ||b0||^2, ||k||^2, b0', k']
% Sliding Window:
% periodically reset the initial measurement y0

%% Data Preparation: y,v,K,dt
K = length(t);
dt = t(2) - t(1);
v = vel;
Iv = zeros(3,K);
y = zeros(1,K);

win_len = 24;

%% Initialization
bias = zeros(3,K);
bias(:,1) = [0.053;-0.011;0.026];

b0 = zeros(3,K);
b0(:,1) = bias(:,1);
k = zeros(3,K);
k(:,1) = [4e-4;1e-4;-3e-4];

x0 = [gtd(1:3,1)', gtd(1:3,1)'*b0(:,1), gtd(1:3,1)'*k(:,1), b0(:,1)'*k(:,1), norm(b0(:,1))^2, norm(k(:,1))^2, b0(:,1)', k(:,1)']';
x_esti = zeros(14,K);
x_pre = zeros(14,K);
x_esti(:,1) = x0;
x_pre(:,1) = x_esti(:,1);

%% x_k = A x_k-1 + B v_k
% A = [eye(3), zeros(3,2), -dt*eye(3); zeros(5,3), eye(5)];
B = [dt*eye(3); zeros(11,3)];

%% Parameter
R = 1e1;
P = diag([1e-2*[1,1,1],1e-2,1e-4,1e-6,1e-4,1e-4,1e-4*[1,1,1],1e-12*[1,1,1]]);
Q = B * diag(1e-3*[0.01,0.01,0.01]) * B';
lambda = 1;

for j = 2:win_len
    A = [eye(3), zeros(3,5), -dt*eye(3), -dt*(j-1)*dt*eye(3); zeros(11,3), eye(11)];
    H = [2*Iv(:,j)',-2*(j-1)*dt,-((j-1)*dt)^2,((j-1)*dt)^3,((j-1)*dt)^2,0.25*((j-1)*dt)^2,zeros(1,6)];
    
    x_pre(:,j) = A * x_esti(:,j-1) + B * v(:,j);
    P = A * P * A'+Q;
    kalman_gain = P * H' * ((H*P*H'+R)^(-1));
    x_esti(:,j) = x_pre(:,j) + kalman_gain * (y(j) - H * x_pre(:,j));
    P = 1/lambda * (eye(14)-kalman_gain*H)*P;
    
    bias(:,j) = x_esti(9:11,j) + (j-1)*dt*x_esti(12:14,j);
end

%% KF_SA_SW
for i = (win_len+1):K
    if rem(i,win_len)==1
        % reset
        start_i = i-1;
        p0 = x_esti(1:3,i-1);
        bias0 = x_esti(9:11,i-1) + dt*win_len*x_esti(12:14,i-1);
        x_esti(9:11,i-1) = bias0;
        y0 = uwb(i-1).^2;
%         x_esti(4,i-1) = p0' * bias0;
%         x_esti(5,i-1) = p0' * x_esti(12:14,i-1);
%         x_esti(6,i-1) = bias0' * x_esti(12:14,i-1);
%         x_esti(7,i-1) = norm(x_esti(9:11,i-1))^2;
%         x_esti(8,i-1) = norm(x_esti(12:14,i-1))^2;
    end
    Iv(1,i) = dt * trapz(v(1,start_i:i));
    Iv(2,i) = dt * trapz(v(2,start_i:i));
    Iv(3,i) = dt * trapz(v(3,start_i:i));
    y(i) = uwb(i).*uwb(i) - y0 + (Iv(1,i)^2 + Iv(2,i)^2 + Iv(3,i)^2);
    
    A = [eye(3), zeros(3,5), -dt*eye(3), -dt*(i-start_i)*dt*eye(3); zeros(11,3), eye(11)];
    H = [2*Iv(:,i)',-2*(i-start_i)*dt,-((i-start_i)*dt)^2,((i-start_i)*dt)^3,((i-start_i)*dt)^2,0.25*((i-start_i)*dt)^2,zeros(1,6)];
    
    x_pre(:,i) = A * x_esti(:,i-1) + B * v(:,i);
    P = A * P * A'+Q;
    kalman_gain = 1/lambda * P * H' * ((H*P*H'+R)^(-1));
    x_esti(:,i) = x_pre(:,i) + kalman_gain * (y(i) - H * x_pre(:,i));
    P = (eye(14)-kalman_gain*H)*P;
    
    bias(:,i) = x_esti(9:11,i) + (i-start_i)*dt*x_esti(12:14,i);
end

figure(12)
plot(t,bias(1,:),'-r',t,bias(2,:),'g-',t,bias(3,:),'b-','linewidth',1);
title('Estimation Bias','FontName','Times New Roman','FontSize',16);
xlabel('Time (s)','FontName','Times New Roman','FontSize',12);
ylabel('Error (m/s)','FontName','Times New Roman','FontSize',12);
legend('X_{bias}','Y_{bias}','Z_{bias}','FontName','Times New Roman','FontSize',8);
