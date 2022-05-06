% function [gtd,t,a,d] = data_processing3()
%% load data
% clc;
close all;
clearvars;
bag_number = 6;

disp(['Bag Number: ',int2str(bag_number)])

uav = importdata(['data2/uav',int2str(bag_number),'.txt']);
vicon = importdata(['data2/vicon',int2str(bag_number),'.txt']);

% world to vicon
qua = uav(2,5:8);
x = qua(1);
y = qua(2);
z = qua(3);
w = qua(4);
R_imu2world = [1-2*y^2-2*z^2, 2*x*y-2*z*w, 2*x*z+2*y*w;
               2*x*y+2*z*w, 1-2*x^2-2*z^2, 2*y*z-2*x*w;
               2*x*z-2*y*w, 2*y*z+2*x*w, 1-2*x^2-2*y^2];
R_world2vicon = R_imu2world';

%% time stamp
dt = 1/30;

% bag1:good; bag2:not that good; bag3:vicon=0; bag456:too long.

% % bag1
% K = 1350;
% start_uav = 923;
% start_vicon = 491;

% % bag2
% K = 1800;
% start_uav = 1052;
% start_vicon = 200;

% % bag4
% K = 2000;
% start_uav = 538;
% start_vicon = 10;

% % bag5
% K = 6000;
% start_uav = 605;
% start_vicon = 10;

% bag6
K = 6000;
start_uav = 866;
start_vicon = 277;

R_all = cell(1,K);

uav = uav(start_uav:(start_uav+K-1),:)';
vicon = vicon(start_vicon:(start_vicon+K-1),:)';

t = 0:dt:K*dt;
t = t(1:K);

% acc = uav(2:4,:);
qua = uav(5:8,:);
acc = uav(2:4,:);
dis = uav(12,:);

a = zeros(3,K);

for i=1:K
    x = qua(1,i);
    y = qua(2,i);
    z = qua(3,i);
    w = qua(4,i);
    R = [1-2*y^2-2*z^2, 2*x*y-2*z*w, 2*x*z+2*y*w;
         2*x*y+2*z*w, 1-2*x^2-2*z^2, 2*y*z-2*x*w;
         2*x*z-2*y*w, 2*y*z+2*x*w, 1-2*x^2-2*y^2];
    
    R_all{1,i} = R * R_world2vicon;
%     acc(3,i) = acc(3,i) - 9.8;
    a(:,i) = R' * acc(:,i);
    a(3,i) = a(3,i) - 9.8;
    a(:,i) = R_world2vicon' * a(:,i);
end
gtd = vicon(2:10,:);

%% S-G filt
A = [eye(3),dt*eye(3);zeros(3,3),eye(3)];
B = [0.5*dt^2*eye(3);dt*eye(3)];

x0 = zeros(6,K);
x = zeros(6,K);
x0 = gtd(1:6,1);
x(:,1) = gtd(1:6,1);

%% Rotation Matrix Check
for i=2:K
    x0(:,i) = A*x0(:,i-1) + B*acc(:,i);
    x(:,i) = A*x(:,i-1) + B*a(:,i);
end
for i=1:3
    figure(i)
    plot(t,x(i,:),'r-',t,gtd(i,:),'b-');
%     plot(t,x0(i,:),'g-',t,x(i,:),'r-',t,gtd(i,:),'b-');
%     plot(t,a(i,:),'g-',t,acc(i,:),'r-',t,gtd(i+6,:),'b-');
end
error0 = a - gtd(7:9,:);
error0_norm = sqrt(mean(error0(1,:).^2 + error0(2,:).^2 + error0(3,:).^2));
disp(error0_norm);
acc = a;

save('bag6.mat' ,'t','dis','acc','gtd');



