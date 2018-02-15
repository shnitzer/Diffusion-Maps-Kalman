function plot_rmse_err( theta1Err, theta2Err, avgSNR, algNames )
%PLOT_RMSE_ERR plots the RMSE of the different algorithms vs. SNR

[SNR,snrSrtI] = sort(mean(avgSNR,2));

figure
subplot(1,2,1); errorbar(mean(SNR,2),mean(theta1Err(snrSrtI,1,:),3),std(theta1Err(snrSrtI,1,:),[],3),'g','LineWidth',1.5); hold on
subplot(1,2,1); errorbar(mean(SNR,2),mean(theta1Err(snrSrtI,2,:),3),std(theta1Err(snrSrtI,2,:),[],3),':m','LineWidth',1);
subplot(1,2,1); errorbar(mean(SNR,2),mean(theta1Err(snrSrtI,3,:),3),std(theta1Err(snrSrtI,3,:),[],3),'-.r','LineWidth',1);
subplot(1,2,1); errorbar(mean(SNR,2),mean(theta1Err(snrSrtI,4,:),3),std(theta1Err(snrSrtI,4,:),[],3),'y','LineWidth',1);
hold off;
set(gca,'xscale','log'); set(gca,'yscale','log'); grid on;
xlabel('SNR (log scale)','FontSize',14); title('$$\theta^{(1)}_n$$','Interpreter','Latex','FontSize',14); ylabel('RMSE (log scale)','FontSize',14)

subplot(1,2,2); errorbar(mean(SNR,2),mean(theta2Err(snrSrtI,1,:),3),std(theta2Err(snrSrtI,1,:),[],3),'g','LineWidth',1.5); hold on
subplot(1,2,2); errorbar(mean(SNR,2),mean(theta2Err(snrSrtI,2,:),3),std(theta2Err(snrSrtI,2,:),[],3),':m','LineWidth',1);
subplot(1,2,2); errorbar(mean(SNR,2),mean(theta2Err(snrSrtI,3,:),3),std(theta2Err(snrSrtI,3,:),[],3),'-.r','LineWidth',1);
subplot(1,2,2); errorbar(mean(SNR,2),mean(theta2Err(snrSrtI,4,:),3),std(theta2Err(snrSrtI,4,:),[],3),'y','LineWidth',1);
hold off;
set(gca,'xscale','log'); set(gca,'yscale','log'); grid on;
xlabel('SNR (log scale)','FontSize',14); title('$$\theta^{(2)}_n$$','Interpreter','Latex','FontSize',14); ylabel('RMSE (log scale)','FontSize',14)

legend(algNames,'FontSize',12);


end

