function [gtd, t, v, d] = data_processing()
%% load data
% clc;
close all;
clearvars;
bag_number = 3;

uav1 = importdata(['data/bag',int2str(bag_number),'_1.txt']);
uav2 = importdata(['data/bag',int2str(bag_number),'_2.txt']);
uav3 = importdata(['data/bag',int2str(bag_number),'_3.txt']);
uav4 = importdata(['data/bag',int2str(bag_number),'_4.txt']);
vicon = importdata(['data/bag',int2str(bag_number),'_vicon.txt']);

%% time stamp
dt = 0.04;

% % bag1
% K = 2700;
% start_1 = 34;
% start_2 = 22;
% start_3 = 43;
% start_4 = 56;
% start_vicon = 135;

% % bag2
% K = 4000;
% start_1 = 105;
% start_2 = 36;
% start_3 = 118;
% start_4 = 119;
% start_vicon = 205; % not set yet

% bag3
K = 2000;
start_1 = 20;
start_2 = 55;
start_3 = 107;
start_4 = 112;
start_vicon = 205; % not set yet

% length=K+1 is to get velocity
uav1 = uav1(start_1:(start_1+K),:);
uav2 = uav2(start_2:(start_2+K),:);
uav3 = uav3(start_3:(start_3+K),:);
uav4 = uav4(start_4:(start_4+K),:);

anchor = uav1;
uav = uav3;

t = 0:dt:K*dt;
t = t(1:K);
v_anchor = (anchor(2:K+1,8:10) - anchor(1:K,8:10))'./dt;
v_uav = (uav(2:K+1,8:10) - uav(1:K,8:10))'./dt;
% v = v_uav - v_anchor;
v = (v_uav - v_anchor) + [0.001;-0.002;0.0005]; % add bias
% depends on the id of anchor and uav
d = uav(1:K,5)';
gtd(1:6,:) = (vicon(start_vicon:start_vicon+K-1,14:19)-vicon(start_vicon:start_vicon+K-1,2:7))';

%% S-G filt
d_raw = d;
order = 1;
framelen = 31;
d = sgolayfilt(d,order,framelen);

v_raw = v;
order = 3;
framelen = 21;
v(1,:) = sgolayfilt(v(1,:),order,framelen);
v(2,:) = sgolayfilt(v(2,:),order,framelen);
v(3,:) = sgolayfilt(v(3,:),order,framelen);

%% test the accurancy of UWB
% d_gtd = sqrt(gtd(1,:).^2 + gtd(2,:).^2 + gtd(3,:).^2);
% figure(6)
% plot(t,d_raw,'g--',t,d_gtd,'b-',t,d,'b--','linewidth',1);
% h1 = legend('uwb_{raw}','gtd','uwb_{filt}','FontName','Times New Roman','FontSize',8,'Location','northwest','NumColumns',1);
% title('SG filt for UWB','FontName','Times New Roman','FontSize',12);
% xlabel('Time','FontName','Times New Roman','FontSize',12);
% ylabel('Distance','FontName','Times New Roman','FontSize',12);
% set(h1,'Orientation','horizon','Box','on');

% d = d_gtd;
% d = d + 0.1; % d have a bias around -0.1m
% plot(t,gtd(4,:),'b-',t,v(1,:),'b--','linewidth',1);

%% test the accurancy of V
% v_gtd = gtd(4:6,:);
% figure(7)
% for i=1:3
%     subplot(3,1,i)
%     plot(t,v_raw(i,:),'g--',t,v_gtd(i,:),'b-',t,v(i,:),'b--','linewidth',1);
%     h1 = legend('V_{raw}','gtd','V_{filt}','FontName','Times New Roman','FontSize',8,'Location','northwest','NumColumns',1);
%     if i == 1
%         title('SG filt for Velocity','FontName','Times New Roman','FontSize',12);
%     end
%     xlabel('Time','FontName','Times New Roman','FontSize',12);
%     ylabel('Distance','FontName','Times New Roman','FontSize',12);
%     set(h1,'Orientation','horizon','Box','on');
% end
