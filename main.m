%% MAIN
% [gtd, t, v, d] = data_processing();
[gtd, t, v, d] = data_processing2();

close all;

% [x_esti,x_predict] = vbakf_r(gtd, t, v, d);
% [error_xyz, error] = result(x_esti, gtd, v, d, t, 4);
% [error_xyz_p, error_p] = result(x_predict, gtd, v, d, t, 5);
% disp(['VBAKF-R ErrorX: ',num2str(error_xyz(1)),'  ErrorY: ',num2str(error_xyz(2)),'  ErrorZ: ',num2str(error_xyz(3)),'  ErrorTotal: ',num2str(error_xyz(4))]);
% disp(['Predict ErrorX: ',num2str(error_xyz_p(1)),'  ErrorY: ',num2str(error_xyz_p(2)),'  ErrorZ: ',num2str(error_xyz_p(3)),'  ErrorTotal: ',num2str(error_xyz_p(4))]);

% [x_kf,x_predict] = kf(gtd, t, v, d);
% [error_xyz_kf, error_kf] = result(x_kf, gtd, v, d, t, 6);
% disp(['KF ErrorX: ',num2str(error_xyz_kf(1)),'  ErrorY: ',num2str(error_xyz_kf(2)),'  ErrorZ: ',num2str(error_xyz_kf(3)),'  ErrorTotal: ',num2str(error_xyz_kf(4))]);

[x_kf_sa] = kf_sa(gtd, t, v, d);
[error_xyz_kf_sa, error_kf_sa] = result(x_kf_sa(1:3,:), gtd, v, d, t, 12);
disp(['KF_SA ErrorX: ',num2str(error_xyz_kf_sa(1)),'  ErrorY: ',num2str(error_xyz_kf_sa(2)),'  ErrorZ: ',num2str(error_xyz_kf_sa(3)),'  ErrorTotal: ',num2str(error_xyz_kf_sa(4))]);

[x_kf_sa3] = kf_sa3(gtd, t, v, d);
[error_xyz_kf_sa3, error_kf_sa3] = result(x_kf_sa3(1:3,:), gtd, v, d, t, 13);
disp(['KF_SA3 ErrorX: ',num2str(error_xyz_kf_sa3(1)),'  ErrorY: ',num2str(error_xyz_kf_sa3(2)),'  ErrorZ: ',num2str(error_xyz_kf_sa3(3)),'  ErrorTotal: ',num2str(error_xyz_kf_sa3(4))]);


% [x_kf_sa_sw] = kf_sa_sw(gtd, t, v, d);
% [error_xyz_kf_sa_sw, error_kf_sa_sw] = result(x_kf_sa_sw(1:3,:), gtd, v, d, t, 14);
% disp(['KF_SA_SW ErrorX: ',num2str(error_xyz_kf_sa_sw(1)),'  ErrorY: ',num2str(error_xyz_kf_sa_sw(2)),'  ErrorZ: ',num2str(error_xyz_kf_sa_sw(3)),'  ErrorTotal: ',num2str(error_xyz_kf_sa_sw(4))]);

% [x_r_sa_sw] = r_sa_sw(gtd, t, v, d);
% [error_xyz_r_sa_sw, error_r_sa_sw] = result(x_r_sa_sw(1:3,:), gtd, v, d, t, 14);
% disp(['R_SA_SW ErrorX: ',num2str(error_xyz_r_sa_sw(1)),'  ErrorY: ',num2str(error_xyz_r_sa_sw(2)),'  ErrorZ: ',num2str(error_xyz_r_sa_sw(3)),'  ErrorTotal: ',num2str(error_xyz_r_sa_sw(4))]);
% 
% [x_q_sa_sw] = q_sa_sw(gtd, t, v, d);
% [error_xyz_q_sa_sw, error_q_sa_sw] = result(x_q_sa_sw(1:3,:), gtd, v, d, t, 14);
% disp(['Q_SA_SW ErrorX: ',num2str(error_xyz_q_sa_sw(1)),'  ErrorY: ',num2str(error_xyz_q_sa_sw(2)),'  ErrorZ: ',num2str(error_xyz_q_sa_sw(3)),'  ErrorTotal: ',num2str(error_xyz_q_sa_sw(4))]);

% [x_q] = vbakf_q(gtd, t, v, d);
% [error_xyz_q, error_q] = result(x_q, gtd, v, d, t, 7);
% disp(['VBAKF-Q ErrorX: ',num2str(error_xyz_q(1)),'  ErrorY: ',num2str(error_xyz_q(2)),'  ErrorZ: ',num2str(error_xyz_q(3)),'  ErrorTotal: ',num2str(error_xyz_q(4))]);

% [x_pm] = vbakf_pm(gtd, t, v, d);3
% [error_xyz_pm, error_pm] = result(x_pm, gtd, v, d, t, 8);
% disp(['VBAKF-PM ErrorX: ',num2str(error_xyz_pm(1)),'  ErrorY: ',num2str(error_xyz_pm(2)),'  ErrorZ: ',num2str(error_xyz_pm(3)),'  ErrorTotal: ',num2str(error_xyz_pm(4))]);

% [x_r_sa] = vbakf_r_sa(gtd, t, v, d);
% [error_xyz_r_sa, error_r_sa] = result(x_r_sa(1:3,:), gtd, v, d, t, 13);
% disp(['KF_R_SA ErrorX: ',num2str(error_xyz_r_sa(1)),'  ErrorY: ',num2str(error_xyz_r_sa(2)),'  ErrorZ: ',num2str(error_xyz_r_sa(3)),'  ErrorTotal: ',num2str(error_xyz_r_sa(4))]);

% [x_kf_sa2] = kf_sa2(gtd, t, v, d);
% [error_xyz_kf_sa2, error_kf_sa2] = result(x_kf_sa2(1:3,:), gtd, v, d, t, 15);
% disp(['KF_SA_2 ErrorX: ',num2str(error_xyz_kf_sa2(1)),'  ErrorY: ',num2str(error_xyz_kf_sa2(2)),'  ErrorZ: ',num2str(error_xyz_kf_sa2(3)),'  ErrorTotal: ',num2str(error_xyz_kf_sa2(4))]);

% [x_kf_sa_sw2] = kf_sa_sw2(gtd, t, v, d);
% [error_xyz_kf_sa_sw2, error_kf_sa_sw2] = result(x_kf_sa_sw2(1:3,:), gtd, v, d, t, 13);
% disp(['KF_SA_SW_2 ErrorX: ',num2str(error_xyz_kf_sa_sw2(1)),'  ErrorY: ',num2str(error_xyz_kf_sa_sw2(2)),'  ErrorZ: ',num2str(error_xyz_kf_sa_sw2(3)),'  ErrorTotal: ',num2str(error_xyz_kf_sa_sw2(4))]);


%% RESULT
% figure(1)
% for i = 1:3
%     subplot(3,1,i)
%     axis = i;
%     plot(t,gtd(axis,:),'b-',t,x_kf_sa(axis,:),'b--',t,x_predict(axis,:),'g-','linewidth',1);
%     hold on
%     h1 = legend('gtd','kf + sa','predict','FontName','Times New Roman','FontSize',12,'Location','northeast','NumColumns',1);
%     xlabel('t (s)','FontName','Times New Roman','FontSize',16);
%     if i ==1
%         ylabel('$r_x$ (m)','FontName','Times New Roman','Interpreter','latex','FontSize',16);
%     end
%     if i ==2
%         ylabel('$r_y$ (m)','FontName','Times New Roman','Interpreter','latex','FontSize',16);
%     end
%     if i ==3
%         ylabel('$r_z$ (m)','FontName','Times New Roman','Interpreter','latex','FontSize',16);
%     end
%     set(h1,'Orientation','horizon','Box','on');
% end
% 
% figure(2)
% plot3(x_kf_sa(1,:),x_kf_sa(2,:),x_kf_sa(3,:),'b--',gtd(1,:),gtd(2,:),gtd(3,:),'b-','linewidth',1);
% set(gca,'linewidth',0.5,'FontSize',8,'FontName','Times New Roman');
% legend('r_{esti}','r_{gtd}','FontName','Times New Roman','FontSize',12);
% ylabel('Trajectory (m)','FontName','Times New Roman','FontSize',14);

% error1 = gtd(1:3,:) - x_kf;
% error_norm1 = sqrt(error1(1,:).^2 + error1(2,:).^2 + error1(3,:).^2);
% error2 = gtd(1:3,:) - x_kf_sa(1:3,:);
% error_norm2 = sqrt(error2(1,:).^2 + error2(2,:).^2 + error2(3,:).^2);
% figure(8)
% plot(t,error_norm1,'b',t,error_norm2,'k','linewidth',1);
% title('RMSE','FontName','Times New Roman','FontSize',16);
% xlabel('Time','FontName','Times New Roman','FontSize',12);
% ylabel('Error','FontName','Times New Roman','FontSize',12);
% legend('KF','KF SA','FontName','Times New Roman','FontSize',8);

%% Bias
% bias1 = x_kf_sa(6:8,:);
% bias2 = x_kf_sa_sw(6:8,:);
% figure(20)
% for i=1:3
%     subplot(3,1,i)
%     plot(t,bias_gtd(i,:),'g',t,bias1(i,:),'b',t,bias2(i,:),'r','linewidth',1);
%     xlabel('Time','FontName','Times New Roman','FontSize',12);
%     ylabel('Bias','FontName','Times New Roman','FontSize',12);
%     legend('gtd','without SW','with SW','FontName','Times New Roman','FontSize',8);
% end
% error1 = sqrt(mean((bias1(1,:)-bias_gtd(1,:)).^2 + (bias1(2,:)-bias_gtd(2,:)).^2 + (bias1(3,:)-bias_gtd(3,:)).^2));
% error2 = sqrt(mean((bias2(1,:)-bias_gtd(1,:)).^2 + (bias2(2,:)-bias_gtd(2,:)).^2 + (bias2(3,:)-bias_gtd(3,:)).^2));
% disp(['Bias_Error ',' without SW: ',num2str(error1),' with SW: ',num2str(error2)]);

% bias_gtd = v - gtd(4:6,:);
% % x_kf_sa(6:8,:) = x_kf_sa(6:8,:)-[0.053;-0.012;0.020];
% bias_gtd = 0.1*bias_gtd + x_kf_sa(6:8,:);
% % x_kf_sa_sw(6:8,:) = x_kf_sa_sw(6:8,:) - [0.053;0.000;0.010];
% order = 3;
% framelen = 31;
% bias_gtd(1,:) = sgolayfilt(bias_gtd(1,:),order,framelen);
% bias_gtd(2,:) = sgolayfilt(bias_gtd(2,:),order,framelen);
% bias_gtd(3,:) = sgolayfilt(bias_gtd(3,:),order,framelen);
% 
% figure(21)
% for i=1:3
%     subplot(3,1,i)
%     plot(t,bias_gtd(i,:),'-b',t,x_kf_sa(i+5,:),'b--','linewidth',1);
%     xlabel('t (s)','FontName','Times New Roman','FontSize',14);
%     ylabel('v_b (m/s)','FontName','Times New Roman','FontSize',14);
%     if i==1
%         legend('$x_{bias}$','$\hat{x}_{bias}$','FontName','Times New Roman','FontSize',12,'Interpreter','latex','NumColumns',1);
%     end
%     if i==2
%         legend('$y_{bias}$','$\hat{y}_{bias}$','FontName','Times New Roman','FontSize',12,'Interpreter','latex','NumColumns',1);
%     end
%     if i==3
%         legend('$z_{bias}$','$\hat{z}_{bias}$','FontName','Times New Roman','FontSize',12,'Interpreter','latex','NumColumns',1);
%     end
% end
