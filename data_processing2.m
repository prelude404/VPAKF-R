% function [gtd, t, v, d] = data_processing2()
%% load data
% clc;
close all;
clearvars;
bag_number = 5;

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
% start_vicon = 491; % not set yet

% % bag2
% K = 1800;
% start_uav = 1052;
% start_vicon = 200;

% % bag3
% K = 1500;
% start_uav = 10;
% start_vicon = 10;

% % bag4
% K = 2000;
% start_uav = 538;
% start_vicon = 10;

% % bag5
% K = 2000;
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
vel = uav(9:11,:);
dis = uav(12,:);

v = zeros(3,K);

for i=1:K
    x = qua(1,i);
    y = qua(2,i);
    z = qua(3,i);
    w = qua(4,i);
    R = [1-2*y^2-2*z^2, 2*x*y-2*z*w, 2*x*z+2*y*w;
         2*x*y+2*z*w, 1-2*x^2-2*z^2, 2*y*z-2*x*w;
         2*x*z-2*y*w, 2*y*z+2*x*w, 1-2*x^2-2*y^2];
    
    R_all{1,i} = R * R_world2vicon;
    
    v(:,i) = (R * R_world2vicon)' * vel(:,i);
end

gtd = vicon(2:7,:);
    
%% S-G filt
% order = 1;
% framelen = 21;
% v(1,:) = sgolayfilt(vel(1,:),order,framelen);
% v(2,:) = sgolayfilt(vel(2,:),order,framelen);
% v(3,:) = sgolayfilt(vel(3,:),order,framelen);
% 
% order = 1;
% framelen = 31;
% d = sgolayfilt(dis,order,framelen);

x0 = zeros(3,K);
x = zeros(3,K);
x0(:,1) = gtd(1:3,1);
x(:,1) = gtd(1:3,1);

A = [1,0,0;
     0,1,0;
     0,0,1];

B = [dt,0,0;
     0,dt,0;
     0,0,dt];
%% Rotation Matrix is necessary
for i=2:K
    x0(:,i) = A*x0(:,i-1) + B*vel(:,i);
    x(:,i) = A*x(:,i-1) + B*v(:,i);
end
for i=1:3
    figure(i)
    plot(t,x(i,:),'r-',t,gtd(i,:),'b-');
%     plot(t,x0(i,:),'g-',t,x(i,:),'r-',t,gtd(i,:),'b-');
%     plot(t,vel(i,:),'g-',t,v(i,:),'r-',t,gtd(i+3,:),'b-');
end
error0 = vel - gtd(4:6,:);
error0_norm = sqrt(mean(error0(1,:).^2 + error0(2,:).^2 + error0(3,:).^2));

% order = 2;
% framelen = 25;
% v(1,:) = sgolayfilt(v(1,:),order,framelen);
% v(2,:) = sgolayfilt(v(2,:),order,framelen);
% v(3,:) = sgolayfilt(v(3,:),order,framelen);
% error = v - gtd(4:6,:);
% error_norm = sqrt(mean(error(1,:).^2 + error(2,:).^2 + error(3,:).^2));
% disp(num2str(error_norm));
% plot(t,v(1,:),'r-',t,gtd(4,:),'b-');

%% Bias exists, bias_mean = [0.0374;-0.0088;0.0063];
% for i=1:3
%     bias(i) = mean(error0(i,:));
%     bias(i+3) = mean(error(i,:));
% end

%% UWB output check
% y_gtd = sqrt(gtd(1,:).^2 + gtd(2,:).^2 + gtd(3,:).^2);
% plot(t,dis,'r-',t,y_gtd,'b-');

% order = 1;
% framelen = 29;
% d = sgolayfilt(dis,order,framelen);
% error_d = sqrt(mean((d-y_gtd).^2));
% disp(num2str(error_d));

