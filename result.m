function [error_xyz,error] = result(x_esti, gtd, vel, uwb, t, fig_num)
figure(fig_num)
set(gcf,'Position',[100,20,600,600]);
subplot(3,1,1)
plot3(x_esti(1,:),x_esti(2,:),x_esti(3,:),'b--',gtd(1,:),gtd(2,:),gtd(3,:),'b-','linewidth',1);
set(gca,'linewidth',0.5,'FontSize',12,'FontName','Times New Roman');
legend('x_{esti}','x_{gtd}','FontName','Times New Roman','FontSize',12);
ylabel('Trajectory (m)','FontName','Times New Roman','FontSize',12);
% title('AKF','FontName','Times New Roman','FontSize',12);
subplot(3,1,2)
plot(t,gtd(1,:),'r-',t,gtd(2,:),'m-',t,gtd(3,:),'b-','linewidth',1);
hold on
plot(t,x_esti(1,:),'r--',t,x_esti(2,:),'m--',t,x_esti(3,:),'b--','linewidth',1);
h1 = legend('x_{gtd}','y_{gtd}','z_{gtd}','x_{esti}','y_{esti}','z_{esti}','Location','northwest','FontName','Times New Roman','FontSize',8);
xlabel('t (s)','FontName','Times New Roman','FontSize',16);
ylabel('r (m)','FontName','Times New Roman','FontSize',16);
set(h1,'Orientation','horizon','Box','on');
set(gca,'linewidth',0.5,'FontSize',12,'FontName','Times New Roman');
subplot(3,1,3)
error = gtd(1:3,:) - x_esti;
error_norm = sqrt(error(1,:).^2 + error(2,:).^2 + error(3,:).^2);
plot(t,error_norm,'k','linewidth',1);
xlabel('t (s)','FontName','Times New Roman','FontSize',16);
ylabel('RMSE (m)','FontName','Times New Roman','FontSize',16);
legend('error','FontName','Times New Roman','FontSize',12);
set(gca,'linewidth',0.5,'FontSize',12,'FontName','Times New Roman');
% hovering part to calculate error
% bag1
% starti = 1250;
% endi = 2250;
% bag2
% starti = 2250;
% endi = 3375;

error = gtd(1:3,:) - x_esti;
error_xyz = zeros(1,4);
error_xyz(1) = sqrt(mean(error(1,:).^2));
error_xyz(2) = sqrt(mean(error(2,:).^2));
error_xyz(3) = sqrt(mean(error(3,:).^2));
error_xyz(4) = sqrt(mean(error(1,:).^2 + error(2,:).^2 + error(3,:).^2));
