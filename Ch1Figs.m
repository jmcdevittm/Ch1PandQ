%Figures derived from data created in PrecipAndQAnalysis
close all

%% Colors and fontsizes
pColor = [31 120 180]./255;
petColor = [51,160,44]./255;
rColor = [227,26,28]./255;
r2Color = [128,177,211]./255;
sColor = [255,127,0]./255;
sDelColor = [253,191,111]./255;
DWColor = [152,78,163]./255;

myFS = 16;

%% Mean Spartanburg precipitation and PET
oneYear = 1:12;

subplot(3,1,1)
hold on
bar(oneYear,SPGmonMeanPrecip,'EdgeColor',pColor,'FaceColor',pColor);
errorbar(oneYear,SPGmonMeanPrecip,SPGmonSDPrecip,'k.')
xlabel('Month','fontsize',myFS)
ylabel('Mean Precip (mm)','fontsize',myFS)
title('Spartanburg Airport 1950 - present','fontsize',myFS)
hold off

subplot(3,1,2)
hold on
bar(oneYear,SPGhistPET,'EdgeColor',petColor,'FaceColor',petColor);
errorbar(oneYear,SPGhistPET,SPGhistPETsd,'k.')
xlabel('Month','fontsize',myFS)
ylabel('Mean PET (mm)','fontsize',myFS)
hold off

subplot(3,1,3)
bar(oneYear,SPGmonMeanPrecip-SPGhistPET,'EdgeColor',r2Color,'FaceColor',r2Color);
xlabel('Month','fontsize',myFS)
ylabel('Mean P-PET (mm)','fontsize',myFS)
ylim([-50 150])
set(gca,'box','off')

%% Monthly Water Balance
figure

myYLim = [0 350];
set(gcf,'defaultAxesColorOrder',[0 0 0; 0 0 0])

subplot(4,1,1)
bar(allDataMonthly.datetime,allDataMonthly.precip,'EdgeColor',pColor,'FaceColor',pColor);
ylabel('Precip (mm)','fontsize',myFS)
ylim(myYLim)

subplot(4,1,2)
bar(allDataMonthly.datetime,allDataMonthly.PET,'EdgeColor',petColor,'FaceColor',petColor)
ylabel('PET (mm)','fontsize',myFS)
ylim(myYLim)

subplot(4,1,3)
yyaxis left
bar(allDataMonthly.datetime,allDataMonthly.runoff,'EdgeColor',rColor,'FaceColor',rColor);
ylabel('Runoff (mm)','fontsize',myFS)
ylim(myYLim)

yyaxis right
scatter(allDataMonthly.datetime,allDataMonthly.RR,'*','MarkerEdgeColor','k')
ylabel('RR','fontsize',myFS)
legend('Runoff','RR')

subplot(4,1,4)
yyaxis right
plot(allDataMonthly.datetime,allDataMonthly.relS,'Color',sColor,'LineWidth',2)
ylabel('Net S (mm)','fontsize',myFS)

yyaxis left
bar(allDataMonthly.datetime,allDataMonthly.delS,'EdgeColor',sDelColor,'FaceColor',sDelColor)
ylabel('\Delta S (mm)','fontsize',myFS)

legend('\DeltaS','Net S')

%% Bivariate Plots

figure
scatter(allDataMonthly.precip,allDataMonthly.runoff,60,'filled','MarkerEdgeColor','none','MarkerFaceColor',pColor)
set(gca,'Fontsize',myFS)
xlabel('Monthly Precip (mm)','fontsize',myFS)
ylabel('Monthly Runoff (mm)','fontsize',myFS)

figure
scatter(allDataMonthly.PET,allDataMonthly.runoff,60,'filled','MarkerEdgeColor','none','MarkerFaceColor',petColor)
set(gca,'Fontsize',myFS)
xlabel('Monthly PET (mm)','fontsize',myFS)
ylabel('Monthly Runoff (mm)','fontsize',myFS)

figure
scatter(allDataMonthly.relS,allDataMonthly.runoff,(100*(allDataMonthly.RR+0.1)),'filled','MarkerEdgeColor','none','MarkerFaceColor',sColor) 
%Added 0.1 to RR because zero values aren't valid for sizing,
%scaled to make points larger
legend('size ~ RR','Location','NW')
xlabel('Monthly NetS (mm)','fontsize',myFS)
ylabel('Monthly Runoff (mm)','fontsize',myFS)
set(gca,'Fontsize',myFS)

%% Daily water balance

figure
set(gcf,'defaultAxesColorOrder',[0 0 0; 0 0 0])
% set(gcf,'Units','inches','Position',[0,0,6.5,7.56])

subplot(6,1,1:2) 

yyaxis left
bar(allDataDaily.datetime,allDataDaily.precip,'EdgeColor',pColor,'FaceColor',pColor); %Precip on first axis
set(gca,'ydir','reverse','ylim',[0 200])
ylabel('Precip (mm/d)','fontsize',myFS)

yyaxis right
plot(allDataDaily.datetime,allDataDaily.runoff,'Color',r2Color,'LineWidth',1) %Runoff on next
set(gca,'ylim',[0 100])
ylabel('Runoff (mm/d)','fontsize',myFS)

DailyXAxisRange = get(gca,'xlim');

%Then plot runoff on logged axis with hyet same as above
subplot(6,1,3:4)

yyaxis left
bar(allDataDaily.datetime,allDataDaily.precip,'EdgeColor',pColor,'FaceColor',pColor); %Precip on first axis
set(gca,'ydir','reverse','ylim',[0 200])
ylabel('Precip (mm/d)','fontsize',myFS)

yyaxis right
semilogy(allDataDaily.datetime,allDataDaily.runoff,'Color',r2Color,'LineWidth',1) %Logged runoff on next
set(gca,'ylim',[0.1 100])
ylabel('Runoff (mm/d)','fontsize',myFS)

%Then water balance residuals
subplot(6,1,5)
plot(allPETMonthly.months(1:36),allPETMonthly.PET(1:36),'Color',petColor,'LineWidth',2) %Exclude right edge
set(gca,'xlim',DailyXAxisRange)
ylabel('PET (mm/d)','fontsize',myFS)

subplot(6,1,6)
plot(allDataDaily.datetime,allDataDaily.relS,'Color',sColor,'LineWidth',2)
set(gca,'xlim',DailyXAxisRange)
ylabel('RelS (mm/d)','fontsize',myFS)

%% WY16 CDF
%Create new timerange
WY2016plusI = timerange('01-Sep-2015 00:00:00','01-Oct-2016 00:00:00');

figure
hold on
plot(allDataDaily.datetime(WY2016plusI),cumsum(allDataDaily.precip(WY2016plusI)),'Color',pColor,'LineWidth',3)
plot(allDataDaily.datetime(WY2016plusI),cumsum(allDataDaily.runoff(WY2016plusI)),'Color',rColor,'LineWidth',3)
plot(allDataDaily.datetime(WY2016plusI),cumsum(allDataDaily.PET(WY2016plusI)),'Color',petColor,'LineWidth',3)
plot(allDataDaily.datetime(WY2016plusI),allDataDaily.relS(WY2016plusI)-min(allDataDaily.relS(WY2016plusI)),'Color',sColor,'LineWidth',3)
legend('Precip','Runoff','PET','relS','location','NW')
set(gca,'Fontsize',myFS)
hold off

%% FDCs

figure
% set(gcf,'Units','inches','Position',[0,0,11,10.5])
hold on
plot(WY2015FDC(:,2),WY2015FDC(:,1),'LineWidth',2)
plot(WY2016FDC(:,2),WY2016FDC(:,1),'LineWidth',2)
plot(WY2017FDC(:,2),WY2017FDC(:,1),'LineWidth',2)
plot(totalFDC(:,2),totalFDC(:,1),'Color','k','LineStyle','--','LineWidth',1)
legend('WY2015','WY2016','WY2017','Mean')
ylim([0.01 100])
set(gca,'YScale','log','Fontsize',myFS)
xlabel('Proportion of year','fontsize',myFS)
ylabel('Runoff (mm/d)','fontsize',myFS)
hold off
