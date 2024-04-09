%Cai
figure;
subplot(2,1,1)
plot(SBfinalt2, SBfinaly2(:,38),'DisplayName','Healthy',"LineWidth",3);
hold on;
plot(SBfinalt_BNd1X2, SBfinaly_BNd1X2(:,38), 'DisplayName','.1X BacNav',"LineWidth",3)
plot(SBfinalt_BNd2X2, SBfinaly_BNd2X2(:,38), 'DisplayName','.2X BacNav',"LineWidth",3)
plot(SBfinalt_BNd4X2, SBfinaly_BNd4X2(:,38), 'DisplayName','.4X BacNav',"LineWidth",3)
plot(SBfinalt_BNd6X2, SBfinaly_BNd6X2(:,38), 'DisplayName','.6X BacNav',"LineWidth",3)
plot(SBfinalt_BNd8X2, SBfinaly_BNd8X2(:,38), 'DisplayName','.8X BacNav',"LineWidth",3)
plot(SBfinalt_BN1X2, SBfinaly_BN1X2(:,38), 'DisplayName','1X BacNav',"LineWidth",3)
title("Healthy Case and BacNav")
xlabel('time (ms)')
ylabel('Ca_i')
legend()
subplot(2,1,2)
plot(SBfinalt2HF, SBfinaly2HF(:,38),'DisplayName','HF',"LineWidth",3);
hold on;
plot(SBfinaltHF_BNd1X2, SBfinalyHF_BNd1X2(:,38), 'DisplayName','HF + .1X BacNav',"LineWidth",3)
plot(SBfinaltHF_BNd2X2, SBfinalyHF_BNd2X2(:,38), 'DisplayName','HF + .2X BacNav',"LineWidth",3)
plot(SBfinaltHF_BNd4X2, SBfinalyHF_BNd4X2(:,38), 'DisplayName','HF + .4X BacNav',"LineWidth",3)
plot(SBfinaltHF_BNd6X2, SBfinalyHF_BNd6X2(:,38), 'DisplayName','HF + .6X BacNav',"LineWidth",3)
plot(SBfinaltHF_BNd8X2, SBfinalyHF_BNd8X2(:,38), 'DisplayName','HF + .8X BacNav',"LineWidth",3)
plot(SBfinaltHF_BN1X2, SBfinalyHF_BN1X2(:,38), 'DisplayName','HF + 1X BacNav',"LineWidth",3)
title("Heart Failure Case and BacNav")
xlabel('time (ms)')
ylabel('Ca_i')
legend()





%CaSR
figure;
subplot(2,1,1)
plot(SBfinalt2, SBfinaly2(:,31),'DisplayName','Healthy',"LineWidth",3);
hold on;
plot(SBfinalt_BNd1X2, SBfinaly_BNd1X2(:,31), 'DisplayName','.1X BacNav',"LineWidth",3)
plot(SBfinalt_BNd2X2, SBfinaly_BNd2X2(:,31), 'DisplayName','.2X BacNav',"LineWidth",3)
plot(SBfinalt_BNd4X2, SBfinaly_BNd4X2(:,31), 'DisplayName','.4X BacNav',"LineWidth",3)
plot(SBfinalt_BNd6X2, SBfinaly_BNd6X2(:,31), 'DisplayName','.6X BacNav',"LineWidth",3)
plot(SBfinalt_BNd8X2, SBfinaly_BNd8X2(:,31), 'DisplayName','.8X BacNav',"LineWidth",3)
plot(SBfinalt_BN1X2, SBfinaly_BN1X2(:,31), 'DisplayName','1X BacNav',"LineWidth",3)
title("Healthy Case and BacNav: CaSR")
xlabel('time (ms)')
ylabel('CaSR')
legend()
subplot(2,1,2)
plot(SBfinalt2HF, SBfinaly2HF(:,31),'DisplayName','HF',"LineWidth",3);
hold on;
plot(SBfinaltHF_BNd1X2, SBfinalyHF_BNd1X2(:,31), 'DisplayName','HF + .1X BacNav',"LineWidth",3)
plot(SBfinaltHF_BNd2X2, SBfinalyHF_BNd2X2(:,31), 'DisplayName','HF + .2X BacNav',"LineWidth",3)
plot(SBfinaltHF_BNd4X2, SBfinalyHF_BNd4X2(:,31), 'DisplayName','HF + .4X BacNav',"LineWidth",3)
plot(SBfinaltHF_BNd6X2, SBfinalyHF_BNd6X2(:,31), 'DisplayName','HF + .6X BacNav',"LineWidth",3)
plot(SBfinaltHF_BNd8X2, SBfinalyHF_BNd8X2(:,31), 'DisplayName','HF + .8X BacNav',"LineWidth",3)
plot(SBfinaltHF_BN1X2, SBfinalyHF_BN1X2(:,31), 'DisplayName','HF + 1X BacNav',"LineWidth",3)
title("Heart Failure Case and BacNav: CaSR")
xlabel('time (ms)')
ylabel('CaSR')
legend()




%I_NCX
figure;
subplot(2,1,1)
plot(SBfinalt2, SBfinalcurrents2(:,6),'DisplayName','Healthy',"LineWidth",3);
hold on;
plot(SBfinalt_BNd1X2, SBfinalcurrents_BNd1X2(:,6), 'DisplayName','.1X BacNav',"LineWidth",3)
plot(SBfinalt_BNd2X2, SBfinalcurrents_BNd2X2(:,6), 'DisplayName','.2X BacNav',"LineWidth",3)
plot(SBfinalt_BNd4X2, SBfinalcurrents_BNd4X2(:,6), 'DisplayName','.4X BacNav',"LineWidth",3)
plot(SBfinalt_BNd6X2, SBfinalcurrents_BNd6X2(:,6), 'DisplayName','.6X BacNav',"LineWidth",3)
plot(SBfinalt_BNd8X2, SBfinalcurrents_BNd8X2(:,6), 'DisplayName','.8X BacNav',"LineWidth",3)
plot(SBfinalt_BN1X2, SBfinalcurrents_BN1X2(:,6), 'DisplayName','1X BacNav',"LineWidth",3)
title("Healthy Case and BacNav: NaCa exchanger")
xlabel('time (ms)')
ylabel('INCX')
legend()
subplot(2,1,2)
plot(SBfinalt2HF, SBfinalcurrents2HF(:,6),'DisplayName','HF',"LineWidth",3);
hold on;
plot(SBfinaltHF_BNd1X2, SBfinalcurrentsHF_BNd1X2(:,6), 'DisplayName','HF + .1X BacNav',"LineWidth",3)
plot(SBfinaltHF_BNd2X2, SBfinalcurrentsHF_BNd2X2(:,6), 'DisplayName','HF + .2X BacNav',"LineWidth",3)
plot(SBfinaltHF_BNd4X2, SBfinalcurrentsHF_BNd4X2(:,6), 'DisplayName','HF + .4X BacNav',"LineWidth",3)
plot(SBfinaltHF_BNd6X2, SBfinalcurrentsHF_BNd6X2(:,6), 'DisplayName','HF + .6X BacNav',"LineWidth",3)
plot(SBfinaltHF_BNd8X2, SBfinalcurrentsHF_BNd8X2(:,6), 'DisplayName','HF + .8X BacNav',"LineWidth",3)
plot(SBfinaltHF_BN1X2, SBfinalcurrentsHF_BN1X2(:,6), 'DisplayName','HF + 1X BacNav',"LineWidth",3)
title("Heart Failure Case and BacNav: NaCa exchanger")
xlabel('time (ms)')
ylabel('INCX')
legend()


% 
% save("TianyuWorkingFinal.mat",'SBfinalt2','SBfinalt_BNd1X2','SBfinalt_BNd2X2','SBfinalt_BNd4X2','SBfinalt_BNd6X2','SBfinalt_BNd8X2','SBfinalt_BN1X2',...
%     'SBfinaly2','SBfinaly_BNd1X2','SBfinaly_BNd2X2','SBfinaly_BNd4X2','SBfinaly_BNd6X2','SBfinaly_BNd8X2','SBfinaly_BN1X2',...
%     'SBfinalcurrents2','SBfinalcurrents_BNd1X2','SBfinalcurrents_BNd2X2','SBfinalcurrents_BNd4X2','SBfinalcurrents_BNd6X2','SBfinalcurrents_BNd8X2','SBfinalcurrents_BN1X2',...
%     'SBfinalt2HF','SBfinaltHF_BNd1X2','SBfinaltHF_BNd2X2','SBfinaltHF_BNd4X2','SBfinaltHF_BNd6X2','SBfinaltHF_BNd8X2','SBfinaltHF_BN1X2',...
%     'SBfinaly2HF','SBfinalyHF_BNd1X2','SBfinalyHF_BNd2X2','SBfinalyHF_BNd4X2','SBfinalyHF_BNd6X2','SBfinalyHF_BNd8X2','SBfinalyHF_BN1X2',...
%     'SBfinalcurrents2HF','SBfinalcurrentsHF_BNd1X2','SBfinalcurrentsHF_BNd2X2','SBfinalcurrentsHF_BNd4X2','SBfinalcurrentsHF_BNd6X2','SBfinalcurrentsHF_BNd8X2','SBfinalcurrentsHF_BN1X2');
% 
% 






%%



%Vm
figure;
subplot(2,1,1)
plot(SBfinalt22, SBfinaly22(:,39),'DisplayName','Healthy',"LineWidth",3);
hold on;
plot(SBfinalt_BNd1X22, SBfinaly_BNd1X22(:,39), 'DisplayName','.1X BacNav',"LineWidth",3)
plot(SBfinalt_BNd2X22, SBfinaly_BNd2X22(:,39), 'DisplayName','.2X BacNav',"LineWidth",3)
plot(SBfinalt_BNd4X22, SBfinaly_BNd4X22(:,39), 'DisplayName','.4X BacNav',"LineWidth",3)
plot(SBfinalt_BNd6X22, SBfinaly_BNd6X22(:,39), 'DisplayName','.6X BacNav',"LineWidth",3)
plot(SBfinalt_BNd8X22, SBfinaly_BNd8X22(:,39), 'DisplayName','.8X BacNav',"LineWidth",3)
plot(SBfinalt_BN1X22, SBfinaly_BN1X22(:,39), 'DisplayName','1X BacNav',"LineWidth",3)
title("Healthy Case and BacNav")
xlabel('time (ms)')
ylabel('Vm')
legend()
subplot(2,1,2)
plot(SBfinalt22HF, SBfinaly22HF(:,39),'DisplayName','HF',"LineWidth",3);
hold on;
plot(SBfinaltHF_BNd1X22, SBfinalyHF_BNd1X22(:,39), 'DisplayName','HF + .1X BacNav',"LineWidth",3)
plot(SBfinaltHF_BNd2X22, SBfinalyHF_BNd2X22(:,39), 'DisplayName','HF + .2X BacNav',"LineWidth",3)
plot(SBfinaltHF_BNd4X22, SBfinalyHF_BNd4X22(:,39), 'DisplayName','HF + .4X BacNav',"LineWidth",3)
plot(SBfinaltHF_BNd6X22, SBfinalyHF_BNd6X22(:,39), 'DisplayName','HF + .6X BacNav',"LineWidth",3)
plot(SBfinaltHF_BNd8X22, SBfinalyHF_BNd8X22(:,39), 'DisplayName','HF + .8X BacNav',"LineWidth",3)
plot(SBfinaltHF_BN1X22, SBfinalyHF_BN1X22(:,39), 'DisplayName','HF + 1X BacNav',"LineWidth",3)
title("Heart Failure Case and BacNav")
xlabel('time (ms)')
ylabel('Vm')
legend()




%Cai
figure;
subplot(2,1,1)
plot(SBfinalt22, SBfinaly22(:,38),'DisplayName','Healthy',"LineWidth",3);
hold on;
plot(SBfinalt_BNd1X22, SBfinaly_BNd1X22(:,38), 'DisplayName','.1X BacNav',"LineWidth",3)
plot(SBfinalt_BNd2X22, SBfinaly_BNd2X22(:,38), 'DisplayName','.2X BacNav',"LineWidth",3)
plot(SBfinalt_BNd4X22, SBfinaly_BNd4X22(:,38), 'DisplayName','.4X BacNav',"LineWidth",3)
plot(SBfinalt_BNd6X22, SBfinaly_BNd6X22(:,38), 'DisplayName','.6X BacNav',"LineWidth",3)
plot(SBfinalt_BNd8X22, SBfinaly_BNd8X22(:,38), 'DisplayName','.8X BacNav',"LineWidth",3)
plot(SBfinalt_BN1X22, SBfinaly_BN1X22(:,38), 'DisplayName','1X BacNav',"LineWidth",3)
title("Healthy Case and BacNav")
xlabel('time (ms)')
ylabel('Ca_i')
legend()
subplot(2,1,2)
plot(SBfinalt22HF, SBfinaly22HF(:,38),'DisplayName','HF',"LineWidth",3);
hold on;
plot(SBfinaltHF_BNd1X22, SBfinalyHF_BNd1X22(:,38), 'DisplayName','HF + .1X BacNav',"LineWidth",3)
plot(SBfinaltHF_BNd2X22, SBfinalyHF_BNd2X22(:,38), 'DisplayName','HF + .2X BacNav',"LineWidth",3)
plot(SBfinaltHF_BNd4X22, SBfinalyHF_BNd4X22(:,38), 'DisplayName','HF + .4X BacNav',"LineWidth",3)
plot(SBfinaltHF_BNd6X22, SBfinalyHF_BNd6X22(:,38), 'DisplayName','HF + .6X BacNav',"LineWidth",3)
plot(SBfinaltHF_BNd8X22, SBfinalyHF_BNd8X22(:,38), 'DisplayName','HF + .8X BacNav',"LineWidth",3)
plot(SBfinaltHF_BN1X22, SBfinalyHF_BN1X22(:,38), 'DisplayName','HF + 1X BacNav',"LineWidth",3)
title("Heart Failure Case and BacNav")
xlabel('time (ms)')
ylabel('Ca_i')
legend()





%CaSR
figure;
subplot(2,1,1)
plot(SBfinalt22, SBfinaly22(:,31),'DisplayName','Healthy',"LineWidth",3);
hold on;
plot(SBfinalt_BNd1X22, SBfinaly_BNd1X22(:,31), 'DisplayName','.1X BacNav',"LineWidth",3)
plot(SBfinalt_BNd2X22, SBfinaly_BNd2X22(:,31), 'DisplayName','.2X BacNav',"LineWidth",3)
plot(SBfinalt_BNd4X22, SBfinaly_BNd4X22(:,31), 'DisplayName','.4X BacNav',"LineWidth",3)
plot(SBfinalt_BNd6X22, SBfinaly_BNd6X22(:,31), 'DisplayName','.6X BacNav',"LineWidth",3)
plot(SBfinalt_BNd8X22, SBfinaly_BNd8X22(:,31), 'DisplayName','.8X BacNav',"LineWidth",3)
plot(SBfinalt_BN1X22, SBfinaly_BN1X22(:,31), 'DisplayName','1X BacNav',"LineWidth",3)
title("Healthy Case and BacNav: CaSR")
xlabel('time (ms)')
ylabel('CaSR')
legend()
subplot(2,1,2)
plot(SBfinalt22HF, SBfinaly22HF(:,31),'DisplayName','HF',"LineWidth",3);
hold on;
plot(SBfinaltHF_BNd1X22, SBfinalyHF_BNd1X22(:,31), 'DisplayName','HF + .1X BacNav',"LineWidth",3)
plot(SBfinaltHF_BNd2X22, SBfinalyHF_BNd2X22(:,31), 'DisplayName','HF + .2X BacNav',"LineWidth",3)
plot(SBfinaltHF_BNd4X22, SBfinalyHF_BNd4X22(:,31), 'DisplayName','HF + .4X BacNav',"LineWidth",3)
plot(SBfinaltHF_BNd6X22, SBfinalyHF_BNd6X22(:,31), 'DisplayName','HF + .6X BacNav',"LineWidth",3)
plot(SBfinaltHF_BNd8X22, SBfinalyHF_BNd8X22(:,31), 'DisplayName','HF + .8X BacNav',"LineWidth",3)
plot(SBfinaltHF_BN1X22, SBfinalyHF_BN1X22(:,31), 'DisplayName','HF + 1X BacNav',"LineWidth",3)
title("Heart Failure Case and BacNav: CaSR")
xlabel('time (ms)')
ylabel('CaSR')
legend()




%I_NCX
figure;
subplot(2,1,1)
plot(SBfinalt22, SBfinalcurrents22(:,6),'DisplayName','Healthy',"LineWidth",3);
hold on;
plot(SBfinalt_BNd1X22, SBfinalcurrents_BNd1X22(:,6), 'DisplayName','.1X BacNav',"LineWidth",3)
plot(SBfinalt_BNd2X22, SBfinalcurrents_BNd2X22(:,6), 'DisplayName','.2X BacNav',"LineWidth",3)
plot(SBfinalt_BNd4X22, SBfinalcurrents_BNd4X22(:,6), 'DisplayName','.4X BacNav',"LineWidth",3)
plot(SBfinalt_BNd6X22, SBfinalcurrents_BNd6X22(:,6), 'DisplayName','.6X BacNav',"LineWidth",3)
plot(SBfinalt_BNd8X22, SBfinalcurrents_BNd8X22(:,6), 'DisplayName','.8X BacNav',"LineWidth",3)
plot(SBfinalt_BN1X22, SBfinalcurrents_BN1X22(:,6), 'DisplayName','1X BacNav',"LineWidth",3)
title("Healthy Case and BacNav: NaCa exchanger")
xlabel('time (ms)')
ylabel('INCX')
legend()
subplot(2,1,2)
plot(SBfinalt22HF, SBfinalcurrents22HF(:,6),'DisplayName','HF',"LineWidth",3);
hold on;
plot(SBfinaltHF_BNd1X22, SBfinalcurrentsHF_BNd1X22(:,6), 'DisplayName','HF + .1X BacNav',"LineWidth",3)
plot(SBfinaltHF_BNd2X22, SBfinalcurrentsHF_BNd2X22(:,6), 'DisplayName','HF + .2X BacNav',"LineWidth",3)
plot(SBfinaltHF_BNd4X22, SBfinalcurrentsHF_BNd4X22(:,6), 'DisplayName','HF + .4X BacNav',"LineWidth",3)
plot(SBfinaltHF_BNd6X22, SBfinalcurrentsHF_BNd6X22(:,6), 'DisplayName','HF + .6X BacNav',"LineWidth",3)
plot(SBfinaltHF_BNd8X22, SBfinalcurrentsHF_BNd8X22(:,6), 'DisplayName','HF + .8X BacNav',"LineWidth",3)
plot(SBfinaltHF_BN1X22, SBfinalcurrentsHF_BN1X22(:,6), 'DisplayName','HF + 1X BacNav',"LineWidth",3)
title("Heart Failure Case and BacNav: NaCa exchanger")
xlabel('time (ms)')
ylabel('INCX')
legend()


%%

% 
% 
% A_Caffeine = { SBfinalt22, SBfinaly22(:,38), SBfinalt22, SBfinaly22(:,31), SBfinalt22, SBfinalcurrents22(:,6), zeros(1000,1),...
%                SBfinalt22HF, SBfinaly22HF(:,38), SBfinalt22HF, SBfinaly22HF(:,31), SBfinalt22HF, SBfinalcurrents22HF(:,6), zeros(1000,1),...
%                SBfinalt_BNd1X22, SBfinaly_BNd1X22(:,38), SBfinalt_BNd1X22, SBfinaly_BNd1X22(:,31), SBfinalt_BNd1X22, SBfinalcurrents_BNd1X22(:,6), zeros(1000,1),...
%                SBfinalt_BNd2X22, SBfinaly_BNd2X22(:,38), SBfinalt_BNd2X22, SBfinaly_BNd2X22(:,31), SBfinalt_BNd2X22, SBfinalcurrents_BNd2X22(:,6), zeros(1000,1),...
%                SBfinalt_BNd4X22, SBfinaly_BNd4X22(:,38), SBfinalt_BNd4X22, SBfinaly_BNd4X22(:,31), SBfinalt_BNd4X22, SBfinalcurrents_BNd4X22(:,6), zeros(1000,1),...
%                SBfinalt_BNd6X22, SBfinaly_BNd6X22(:,38), SBfinalt_BNd6X22, SBfinaly_BNd6X22(:,31), SBfinalt_BNd6X22, SBfinalcurrents_BNd6X22(:,6), zeros(1000,1),...
%                SBfinalt_BNd8X22, SBfinaly_BNd8X22(:,38), SBfinalt_BNd8X22, SBfinaly_BNd8X22(:,31), SBfinalt_BNd8X22, SBfinalcurrents_BNd8X22(:,6), zeros(1000,1),...
%                SBfinalt_BN1X22, SBfinaly_BN1X22(:,38), SBfinalt_BN1X22, SBfinaly_BN1X22(:,31), SBfinalt_BN1X22, SBfinalcurrents_BN1X22(:,6), zeros(1000,1),...
%                SBfinaltHF_BNd1X22, SBfinalyHF_BNd1X22(:,38), SBfinaltHF_BNd1X22, SBfinalyHF_BNd1X22(:,31), SBfinaltHF_BNd1X22, SBfinalcurrentsHF_BNd1X22(:,6), zeros(1000,1),...
%                SBfinaltHF_BNd2X22, SBfinalyHF_BNd2X22(:,38), SBfinaltHF_BNd2X22, SBfinalyHF_BNd2X22(:,31), SBfinaltHF_BNd2X22, SBfinalcurrentsHF_BNd2X22(:,6), zeros(1000,1),...
%                SBfinaltHF_BNd4X22, SBfinalyHF_BNd4X22(:,38), SBfinaltHF_BNd4X22, SBfinalyHF_BNd4X22(:,31), SBfinaltHF_BNd4X22, SBfinalcurrentsHF_BNd4X22(:,6), zeros(1000,1),...
%                SBfinaltHF_BNd6X22, SBfinalyHF_BNd6X22(:,38), SBfinaltHF_BNd6X22, SBfinalyHF_BNd6X22(:,31), SBfinaltHF_BNd6X22, SBfinalcurrentsHF_BNd6X22(:,6), zeros(1000,1),...
%                SBfinaltHF_BNd8X22, SBfinalyHF_BNd8X22(:,38), SBfinaltHF_BNd8X22, SBfinalyHF_BNd8X22(:,31), SBfinaltHF_BNd8X22, SBfinalcurrentsHF_BNd8X22(:,6), zeros(1000,1),...
%                SBfinaltHF_BN1X22, SBfinalyHF_BN1X22(:,38), SBfinaltHF_BN1X22, SBfinalyHF_BN1X22(:,31), SBfinaltHF_BN1X22, SBfinalcurrentsHF_BN1X22(:,6)};

filename = '../SimResults/DADs_VmAndCai.xlsx';
alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ';


dataToExportBACNAV = { SBfinaltKir1, SBfinalyKir1(:,39), timeImportKir1, yImportKir1(:,39), SBfinaltKir1, SBfinalyKir1(:,38), timeImportKir1, yImportKir1(:,38), zeros(1000,1),...
                 SBfinalt_BNd1XKir1, SBfinaly_BNd1XKir1(:,39), timeImport_BNd1XKir1, yImport_BNd1XKir1(:,39), SBfinalt_BNd1XKir1, SBfinaly_BNd1XKir1(:,38), timeImport_BNd1XKir1, yImport_BNd1XKir1(:,38), zeros(1000,1),...
                 SBfinalt_BNd2XKir1, SBfinaly_BNd2XKir1(:,39), timeImport_BNd2XKir1, yImport_BNd2XKir1(:,39), SBfinalt_BNd2XKir1, SBfinaly_BNd2XKir1(:,38), timeImport_BNd2XKir1, yImport_BNd2XKir1(:,38), zeros(1000,1),...
                 SBfinalt_BNd4XKir1, SBfinaly_BNd4XKir1(:,39), timeImport_BNd4XKir1, yImport_BNd4XKir1(:,39), SBfinalt_BNd4XKir1, SBfinaly_BNd4XKir1(:,38), timeImport_BNd4XKir1, yImport_BNd4XKir1(:,38), zeros(1000,1),...
                 SBfinalt_BNd6XKir1, SBfinaly_BNd6XKir1(:,39), timeImport_BNd6XKir1, yImport_BNd6XKir1(:,39), SBfinalt_BNd6XKir1, SBfinaly_BNd6XKir1(:,38), timeImport_BNd6XKir1, yImport_BNd6XKir1(:,38), zeros(1000,1),...
                 SBfinalt_BNd8XKir1, SBfinaly_BNd8XKir1(:,39), timeImport_BNd8XKir1, yImport_BNd8XKir1(:,39), SBfinalt_BNd8XKir1, SBfinaly_BNd8XKir1(:,38), timeImport_BNd8XKir1, yImport_BNd8XKir1(:,38), zeros(1000,1),...
                 SBfinalt_BN1XKir1, SBfinaly_BN1XKir1(:,39), timeImport_BN1XKir1, yImport_BN1XKir1(:,39), SBfinalt_BN1XKir1, SBfinaly_BN1XKir1(:,38), timeImport_BN1XKir1, yImport_BN1XKir1(:,38), zeros(1000,1) };



% Write the arrays to the Excel file, each array in a separate column
for i = 1:length(dataToExportBACNAV)

    if i > 78
        range = [ 'C' alphabet(i-78) num2str(1) ':C' alphabet(i-78) num2str(length(dataToExportBACNAV{i})) ];
    elseif i > 52
        range = [ 'B' alphabet(i-52) num2str(1) ':B' alphabet(i-52) num2str(length(dataToExportBACNAV{i})) ];
    elseif i > 26
        range = [ 'A' alphabet(i-26) num2str(1) ':A' alphabet(i-26) num2str(length(dataToExportBACNAV{i})) ];
    else
        range = [ alphabet(i) num2str(1) ':' alphabet(i) num2str(length(dataToExportBACNAV{i})) ]; 
    end
        

    writematrix(dataToExportBACNAV{i}, filename, 'Range', range );
    % 'Range' specifies the column (A, B, C, etc.) where data is written
    % 'WriteMode' is set to 'append' to add columns for each array
end





%%


tableNCXInhibition = { SBfinaltNCXMult_0, SBfinalyNCXMult_0(:,39), SBfinalcurrentsNCXMult_0(:,6), zeros(1000,1),...
                       SBfinalt_BNd1XNCXMult_0, SBfinaly_BNd1XNCXMult_0(:,39), SBfinalcurrents_BNd1XNCXMult_0(:,6), zeros(1000,1),...
                       SBfinalt_BNd2XNCXMult_0, SBfinaly_BNd2XNCXMult_0(:,39), SBfinalcurrents_BNd2XNCXMult_0(:,6), zeros(1000,1),...
                       SBfinalt_BNd4XNCXMult_0, SBfinaly_BNd4XNCXMult_0(:,39), SBfinalcurrents_BNd4XNCXMult_0(:,6), zeros(1000,1),...
                       SBfinalt_BNd6XNCXMult_0, SBfinaly_BNd6XNCXMult_0(:,39), SBfinalcurrents_BNd6XNCXMult_0(:,6), zeros(1000,1),...
                       SBfinalt_BNd8XNCXMult_0, SBfinaly_BNd8XNCXMult_0(:,39), SBfinalcurrents_BNd8XNCXMult_0(:,6), zeros(1000,1),...
                       SBfinalt_BN1XNCXMult_0, SBfinaly_BN1XNCXMult_0(:,39), SBfinalcurrents_BN1XNCXMult_0(:,6), zeros(1000,1),...
                        SBfinaltNCXMult_d01, SBfinalyNCXMult_d01(:,39), SBfinalcurrentsNCXMult_d01(:,6), zeros(1000,1),...
                       SBfinalt_BNd1XNCXMult_d01, SBfinaly_BNd1XNCXMult_d01(:,39), SBfinalcurrents_BNd1XNCXMult_d01(:,6), zeros(1000,1),...
                       SBfinalt_BNd2XNCXMult_d01, SBfinaly_BNd2XNCXMult_d01(:,39), SBfinalcurrents_BNd2XNCXMult_d01(:,6), zeros(1000,1),...
                       SBfinalt_BNd4XNCXMult_d01, SBfinaly_BNd4XNCXMult_d01(:,39), SBfinalcurrents_BNd4XNCXMult_d01(:,6), zeros(1000,1),...
                       SBfinalt_BNd6XNCXMult_d01, SBfinaly_BNd6XNCXMult_d01(:,39), SBfinalcurrents_BNd6XNCXMult_d01(:,6), zeros(1000,1),...
                       SBfinalt_BNd8XNCXMult_d01, SBfinaly_BNd8XNCXMult_d01(:,39), SBfinalcurrents_BNd8XNCXMult_d01(:,6), zeros(1000,1),...
                       SBfinalt_BN1XNCXMult_d01, SBfinaly_BN1XNCXMult_d01(:,39), SBfinalcurrents_BN1XNCXMult_d01(:,6), zeros(1000,1),...
                        SBfinaltNCXMult_d05, SBfinalyNCXMult_d05(:,39), SBfinalcurrentsNCXMult_d05(:,6), zeros(1000,1),...
                       SBfinalt_BNd1XNCXMult_d05, SBfinaly_BNd1XNCXMult_d05(:,39), SBfinalcurrents_BNd1XNCXMult_d05(:,6), zeros(1000,1),...
                       SBfinalt_BNd2XNCXMult_d05, SBfinaly_BNd2XNCXMult_d05(:,39), SBfinalcurrents_BNd2XNCXMult_d05(:,6), zeros(1000,1),...
                       SBfinalt_BNd4XNCXMult_d05, SBfinaly_BNd4XNCXMult_d05(:,39), SBfinalcurrents_BNd4XNCXMult_d05(:,6), zeros(1000,1),...
                       SBfinalt_BNd6XNCXMult_d05, SBfinaly_BNd6XNCXMult_d05(:,39), SBfinalcurrents_BNd6XNCXMult_d05(:,6), zeros(1000,1),...
                       SBfinalt_BNd8XNCXMult_d05, SBfinaly_BNd8XNCXMult_d05(:,39), SBfinalcurrents_BNd8XNCXMult_d05(:,6), zeros(1000,1),...
                       SBfinalt_BN1XNCXMult_d05, SBfinaly_BN1XNCXMult_d05(:,39), SBfinalcurrents_BN1XNCXMult_d05(:,6), zeros(1000,1),...
                        SBfinaltNCXMult_d1, SBfinalyNCXMult_d1(:,39), SBfinalcurrentsNCXMult_d1(:,6), zeros(1000,1),...
                       SBfinalt_BNd1XNCXMult_d1, SBfinaly_BNd1XNCXMult_d1(:,39), SBfinalcurrents_BNd1XNCXMult_d1(:,6), zeros(1000,1),...
                       SBfinalt_BNd2XNCXMult_d1, SBfinaly_BNd2XNCXMult_d1(:,39), SBfinalcurrents_BNd2XNCXMult_d1(:,6), zeros(1000,1),...
                       SBfinalt_BNd4XNCXMult_d1, SBfinaly_BNd4XNCXMult_d1(:,39), SBfinalcurrents_BNd4XNCXMult_d1(:,6), zeros(1000,1),...
                       SBfinalt_BNd6XNCXMult_d1, SBfinaly_BNd6XNCXMult_d1(:,39), SBfinalcurrents_BNd6XNCXMult_d1(:,6), zeros(1000,1),...
                       SBfinalt_BNd8XNCXMult_d1, SBfinaly_BNd8XNCXMult_d1(:,39), SBfinalcurrents_BNd8XNCXMult_d1(:,6), zeros(1000,1),...
                       SBfinalt_BN1XNCXMult_d1, SBfinaly_BN1XNCXMult_d1(:,39), SBfinalcurrents_BN1XNCXMult_d1(:,6), zeros(1000,1),...
                        SBfinaltNCXMult_d2, SBfinalyNCXMult_d2(:,39), SBfinalcurrentsNCXMult_d2(:,6), zeros(1000,1),...
                       SBfinalt_BNd1XNCXMult_d2, SBfinaly_BNd1XNCXMult_d2(:,39), SBfinalcurrents_BNd1XNCXMult_d2(:,6), zeros(1000,1),...
                       SBfinalt_BNd2XNCXMult_d2, SBfinaly_BNd2XNCXMult_d2(:,39), SBfinalcurrents_BNd2XNCXMult_d2(:,6), zeros(1000,1),...
                       SBfinalt_BNd4XNCXMult_d2, SBfinaly_BNd4XNCXMult_d2(:,39), SBfinalcurrents_BNd4XNCXMult_d2(:,6), zeros(1000,1),...
                       SBfinalt_BNd6XNCXMult_d2, SBfinaly_BNd6XNCXMult_d2(:,39), SBfinalcurrents_BNd6XNCXMult_d2(:,6), zeros(1000,1),...
                       SBfinalt_BNd8XNCXMult_d2, SBfinaly_BNd8XNCXMult_d2(:,39), SBfinalcurrents_BNd8XNCXMult_d2(:,6), zeros(1000,1),...
                       SBfinalt_BN1XNCXMult_d2, SBfinaly_BN1XNCXMult_d2(:,39), SBfinalcurrents_BN1XNCXMult_d2(:,6), zeros(1000,1),...
                        SBfinaltNCXMult_d4, SBfinalyNCXMult_d4(:,39), SBfinalcurrentsNCXMult_d4(:,6), zeros(1000,1),...
                       SBfinalt_BNd1XNCXMult_d4, SBfinaly_BNd1XNCXMult_d4(:,39), SBfinalcurrents_BNd1XNCXMult_d4(:,6), zeros(1000,1),...
                       SBfinalt_BNd2XNCXMult_d4, SBfinaly_BNd2XNCXMult_d4(:,39), SBfinalcurrents_BNd2XNCXMult_d4(:,6), zeros(1000,1),...
                       SBfinalt_BNd4XNCXMult_d4, SBfinaly_BNd4XNCXMult_d4(:,39), SBfinalcurrents_BNd4XNCXMult_d4(:,6), zeros(1000,1),...
                       SBfinalt_BNd6XNCXMult_d4, SBfinaly_BNd6XNCXMult_d4(:,39), SBfinalcurrents_BNd6XNCXMult_d4(:,6), zeros(1000,1),...
                       SBfinalt_BNd8XNCXMult_d4, SBfinaly_BNd8XNCXMult_d4(:,39), SBfinalcurrents_BNd8XNCXMult_d4(:,6), zeros(1000,1),...
                       SBfinalt_BN1XNCXMult_d4, SBfinaly_BN1XNCXMult_d4(:,39), SBfinalcurrents_BN1XNCXMult_d4(:,6), zeros(1000,1),...
                        SBfinaltNCXMult_d5, SBfinalyNCXMult_d5(:,39), SBfinalcurrentsNCXMult_d5(:,6), zeros(1000,1),...
                       SBfinalt_BNd1XNCXMult_d5, SBfinaly_BNd1XNCXMult_d5(:,39), SBfinalcurrents_BNd1XNCXMult_d5(:,6), zeros(1000,1),...
                       SBfinalt_BNd2XNCXMult_d5, SBfinaly_BNd2XNCXMult_d5(:,39), SBfinalcurrents_BNd2XNCXMult_d5(:,6), zeros(1000,1),...
                       SBfinalt_BNd4XNCXMult_d5, SBfinaly_BNd4XNCXMult_d5(:,39), SBfinalcurrents_BNd4XNCXMult_d5(:,6), zeros(1000,1),...
                       SBfinalt_BNd6XNCXMult_d5, SBfinaly_BNd6XNCXMult_d5(:,39), SBfinalcurrents_BNd6XNCXMult_d5(:,6), zeros(1000,1),...
                       SBfinalt_BNd8XNCXMult_d5, SBfinaly_BNd8XNCXMult_d5(:,39), SBfinalcurrents_BNd8XNCXMult_d5(:,6), zeros(1000,1),...
                       SBfinalt_BN1XNCXMult_d5, SBfinaly_BN1XNCXMult_d5(:,39), SBfinalcurrents_BN1XNCXMult_d5(:,6), zeros(1000,1),...
                        SBfinaltNCXMult_d6, SBfinalyNCXMult_d6(:,39), SBfinalcurrentsNCXMult_d6(:,6), zeros(1000,1),...
                       SBfinalt_BNd1XNCXMult_d6, SBfinaly_BNd1XNCXMult_d6(:,39), SBfinalcurrents_BNd1XNCXMult_d6(:,6), zeros(1000,1),...
                       SBfinalt_BNd2XNCXMult_d6, SBfinaly_BNd2XNCXMult_d6(:,39), SBfinalcurrents_BNd2XNCXMult_d6(:,6), zeros(1000,1),...
                       SBfinalt_BNd4XNCXMult_d6, SBfinaly_BNd4XNCXMult_d6(:,39), SBfinalcurrents_BNd4XNCXMult_d6(:,6), zeros(1000,1),...
                       SBfinalt_BNd6XNCXMult_d6, SBfinaly_BNd6XNCXMult_d6(:,39), SBfinalcurrents_BNd6XNCXMult_d6(:,6), zeros(1000,1),...
                       SBfinalt_BNd8XNCXMult_d6, SBfinaly_BNd8XNCXMult_d6(:,39), SBfinalcurrents_BNd8XNCXMult_d6(:,6), zeros(1000,1),...
                       SBfinalt_BN1XNCXMult_d6, SBfinaly_BN1XNCXMult_d6(:,39), SBfinalcurrents_BN1XNCXMult_d6(:,6), zeros(1000,1),...
                        SBfinaltNCXMult_d8, SBfinalyNCXMult_d8(:,39), SBfinalcurrentsNCXMult_d8(:,6), zeros(1000,1),...
                       SBfinalt_BNd1XNCXMult_d8, SBfinaly_BNd1XNCXMult_d8(:,39), SBfinalcurrents_BNd1XNCXMult_d8(:,6), zeros(1000,1),...
                       SBfinalt_BNd2XNCXMult_d8, SBfinaly_BNd2XNCXMult_d8(:,39), SBfinalcurrents_BNd2XNCXMult_d8(:,6), zeros(1000,1),...
                       SBfinalt_BNd4XNCXMult_d8, SBfinaly_BNd4XNCXMult_d8(:,39), SBfinalcurrents_BNd4XNCXMult_d8(:,6), zeros(1000,1),...
                       SBfinalt_BNd6XNCXMult_d8, SBfinaly_BNd6XNCXMult_d8(:,39), SBfinalcurrents_BNd6XNCXMult_d8(:,6), zeros(1000,1),...
                       SBfinalt_BNd8XNCXMult_d8, SBfinaly_BNd8XNCXMult_d8(:,39), SBfinalcurrents_BNd8XNCXMult_d8(:,6), zeros(1000,1),...
                       SBfinalt_BN1XNCXMult_d8, SBfinaly_BN1XNCXMult_d8(:,39), SBfinalcurrents_BN1XNCXMult_d8(:,6), zeros(1000,1),...
                        SBfinaltNCXMult_1, SBfinalyNCXMult_1(:,39), SBfinalcurrentsNCXMult_1(:,6), zeros(1000,1),...
                       SBfinalt_BNd1XNCXMult_1, SBfinaly_BNd1XNCXMult_1(:,39), SBfinalcurrents_BNd1XNCXMult_1(:,6), zeros(1000,1),...
                       SBfinalt_BNd2XNCXMult_1, SBfinaly_BNd2XNCXMult_1(:,39), SBfinalcurrents_BNd2XNCXMult_1(:,6), zeros(1000,1),...
                       SBfinalt_BNd4XNCXMult_1, SBfinaly_BNd4XNCXMult_1(:,39), SBfinalcurrents_BNd4XNCXMult_1(:,6), zeros(1000,1),...
                       SBfinalt_BNd6XNCXMult_1, SBfinaly_BNd6XNCXMult_1(:,39), SBfinalcurrents_BNd6XNCXMult_1(:,6), zeros(1000,1),...
                       SBfinalt_BNd8XNCXMult_1, SBfinaly_BNd8XNCXMult_1(:,39), SBfinalcurrents_BNd8XNCXMult_1(:,6), zeros(1000,1),...
                       SBfinalt_BN1XNCXMult_1, SBfinaly_BN1XNCXMult_1(:,39), SBfinalcurrents_BN1XNCXMult_1(:,6) }; 


createExcelFromCell(tableNCXInhibition, "../SimResults/INCX_inhibition_data.xlsx")





















