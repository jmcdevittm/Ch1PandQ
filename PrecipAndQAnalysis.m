%% Runoff Ratio Analyses of precip and Q for all 3 water years 2015-2017
clearvars -except CalhounData
close all

%Import .csv's of precip for WYs 2015 through 2017 from various folders,
%and then save as .mat to facilitate analysis
% WY2015precip = table2timetable(readtable('CalhounPrecipWY2015.csv'));
% WY2016precip = table2timetable(readtable('CalhounPrecipWY2016.csv'));
% WY2017precip = table2timetable(readtable('WY2017precip.csv'));
% save('AllPrecipData.mat','WY*') %Save as .mat

%% Create Water year Indicess
WY2015I = timerange('01-Oct-2014 00:00:00','01-Oct-2015 00:00:00');
WY2016I = timerange('01-Oct-2015 00:00:00','01-Oct-2016 00:00:00');
WY2017I = timerange('01-Oct-2016 00:00:00','01-Oct-2017 00:00:00');

%% Load Precip and Q/Runoff, standardize to 5 min resolution, and synchronize into same timetable
%Precip from AllPrecipData .mat file, Q/Runoff from CalhounData
load AllPrecipData.mat

%Concatenate to form total record
allPrecip = vertcat(WY2015precip,WY2016precip,WY2017precip);

%Standardize times to 5 min resolution. Do not interpolate Precip
[~,ind] = unique(allPrecip.datetime,'last'); %Find and remove duplicate values
dupI = setdiff(1:height(allPrecip),ind);
allPrecip(dupI,:) = [];
%Put on regular time series. No interp, add nearest value on original ts to new ts
allPrecip = retime(allPrecip,'regular','nearest','TimeStep',minutes(5));

%Q/Runoff
allI = CalhounData(37).data.datetime >= datenum('100114 00:00:00','mmddyy HH:MM:SS') & ...
    CalhounData(37).data.datetime < datenum('100117 00:00:00','mmddyy HH:MM:SS');

allQRunoff = CalhounData(37).data(allI,[1 6 7]);
allQRunoff.datetime = datetime(allQRunoff.datetime,'ConvertFrom','datenum');
allQRunoff = table2timetable(allQRunoff);

%Standardize times to 5 min resolution. Interpolate Q/R
[~,ind] = unique(allQRunoff.datetime); %Find and remove duplicate values
dupI = setdiff(1:height(allQRunoff),ind);
allQRunoff(dupI,:) = [];
allQRunoff = retime(allQRunoff,'regular','linear','TimeStep',minutes(5)); %Put time series on regular time interval

allData = synchronize(allQRunoff,allPrecip); %Synchronize onto same time series

%% Load Deep Well Data, CalhounData(28) and trim to three water years, interp or add noData values
allGWData = CalhounData(28).data(:,1:2);

%convert to timetable
allGWData = timetable(datetime(allGWData.datetime,'ConvertFrom','datenum'),allGWData.level);
allGWData.Properties.VariableNames = {'DWlevel'};

%Put data on regular, 20 min time interval
allGWData = retime(allGWData,'regular','linear','TimeStep',minutes(20));

%Data gaps (Determined manually from plot):
%13-Apr-2016 21:20:00 -> 30-Apr-2016 19:20:00
%05-Aug-2016 05:00:00 -> 10-Aug-2016 13:00:00
%05-Oct-2016 02:20:00 -> 18-Oct-2016 13:20:00
dataGaps =... 
    {timerange('13-Apr-2016 21:20:00','30-Apr-2016 19:20:00','open')...
    timerange('05-Aug-2016 05:00:00','10-Aug-2016 13:00:00','open')...
    timerange('05-Oct-2016 02:20:00','18-Oct-2016 13:20:00','open')};

%Trim to just three water years
allGWData = allGWData(timerange('01-Oct-2014 00:00:00','01-Oct-2017 00:00:00'),:);

%% Calculate Thornwaithe PET in monthly and daily time steps
load meanMonthlyTemp
load dayLength

daylength = daylength * 24;
monthDays = [30 31 30 31 31 28 31 30 31 30 31 31]';

%TW PET estimate 
WY2015HeatIndex = sum((meanMonthTemp(1:12)./5).^1.514);
WY2015alpha = 6.75e-7*WY2015HeatIndex^3 - 7.71e-5*WY2015HeatIndex^2 + 1.792e-2*WY2015HeatIndex + 0.49239;
WY2015PET = 16.*(daylength./12).*(monthDays./30).*(10.*meanMonthTemp(1:12)./WY2015HeatIndex).^WY2015alpha; %mm/month
WY2015PET = WY2015PET./monthDays; %mm/day

WY2016HeatIndex = sum((meanMonthTemp(13:24)./5).^1.514);
WY2016alpha = 6.75e-7*WY2016HeatIndex^3 - 7.71e-5*WY2016HeatIndex^2 + 1.792e-2*WY2016HeatIndex + 0.49239;
WY2016PET = 16.*(daylength./12).*(monthDays./30).*(10.*meanMonthTemp(13:24)./WY2016HeatIndex).^WY2016alpha; %mm/month
WY2016PET = WY2016PET./monthDays; %mm/day

WY2017HeatIndex = sum((meanMonthTemp(25:36)./5).^1.514);
WY2017alpha = 6.75e-7*WY2017HeatIndex^3 - 7.71e-5*WY2017HeatIndex^2 + 1.792e-2*WY2017HeatIndex + 0.49239;
WY2017PET = 16.*(daylength./12).*(monthDays./30).*(10.*meanMonthTemp(25:36)./WY2017HeatIndex).^WY2017alpha; %mm/month
WY2017PET = WY2017PET./monthDays; %mm/day

PET = [WY2015PET; WY2016PET; WY2017PET; 0]; %Add zero to macth w/ Oct 17 later when interpolating to daily
months = (datetime(2014,10,1):calmonths(1):datetime(2017,10,1))'; %Create monthly time series
allPETMonthly = timetable(months,PET);
allPETDaily = retime(allPETMonthly,'daily','previous');
allPETDaily(end,:) = []; %Remove row for Oct 1 2017

%% Create HOURLY runoff/q, precip
allDataHourly = allData;
allDataHourly.runoff = allDataHourly.runoff/12; %Convert from mm/hr to just depth of water over 5 min interval
allDataHourly.discharge = allDataHourly.discharge*300; %Convert from L/s to just volume of water over 5 min interval
allDataHourly = retime(allDataHourly,'hourly','sum');

%Create categorical variable for WY (i.e, 15,16,17)
allDataHourly.WY = zeros(height(allDataHourly),1);
allDataHourly.WY(WY2015I) = 15;
allDataHourly.WY(WY2016I) = 16;
allDataHourly.WY(WY2017I) = 17;

%Create categorical variable for growing/dormant seasons (1 = growing; 0 =
%dormant)
allDataHourly.seas = zeros(height(allDataHourly),1);
allDataHourly.seas = month(allDataHourly.datetime);
allDataHourly.seas = double(allDataHourly.seas > 3 & allDataHourly.seas < 11);

%% Create DAILY runoff, precip, etc
allDataDaily = allData;
allDataDaily.runoff = allDataDaily.runoff/12; %Convert from mm/hr to just depth of water over 5 min interval
allDataDaily.discharge = allDataDaily.discharge*300; %Convert from L/s to just volume of water over 5 min interval
allDataDaily = retime(allDataDaily,'daily','sum');

allDataDaily = synchronize(allDataDaily,retime(allGWData,'daily','linear'),'intersection','linear'); %Add GW data to daily data

allDataDaily = synchronize(allDataDaily,allPETDaily); %Add PET to daily data

allDataDaily.delS = zeros(height(allDataDaily),1); %Calculate daily change in storage
allDataDaily.delS = allDataDaily.precip - allDataDaily.runoff - allDataDaily.PET;

allDataDaily.relS = zeros(height(allDataDaily),1); %Calculate relative storage from change in storage
allDataDaily.relS = cumsum(allDataDaily.delS);

%Create categorical variable for WY (i.e, 15,16,17)
allDataDaily.WY = zeros(height(allDataDaily),1);
allDataDaily.WY(WY2015I) = 15;
allDataDaily.WY(WY2016I) = 16;
allDataDaily.WY(WY2017I) = 17;

%Create categorical variable for growing/dormant seasons (1 = growing; 0 =
%dormant)
allDataDaily.seas = zeros(height(allDataDaily),1);
allDataDaily.seas = month(allDataDaily.datetime);
allDataDaily.seas = double(allDataDaily.seas > 3 & allDataDaily.seas < 11);

%% Create MONTHLY runoff, precip, etc
allDataMonthly = allDataDaily(:,1:7); %Exclude categorical values from above
allDataMonthly = retime(allDataMonthly,'monthly','sum');
allDataMonthly.relS = cumsum(allDataMonthly.delS);
allDataMonthly.RR = allDataMonthly.runoff./allDataMonthly.precip;
allDataMonthly.WI = allDataMonthly.precip./allDataMonthly.PET;

%Add GW Data
allDataMonthly.DWlevel = [];
GWtemp = retime(allGWData,'monthly','linear');
allDataMonthly.DWlevel = GWtemp.DWlevel(1:end-1); %Add GW data to daily data

%Create categorical variable for WY (i.e, 15,16,17)
allDataMonthly.WY = zeros(height(allDataMonthly),1);
allDataMonthly.WY(WY2015I) = 15;
allDataMonthly.WY(WY2016I) = 16;
allDataMonthly.WY(WY2017I) = 17;

%Create categorical variable for growing/dormant seasons (1 = growing; 0 =
%dormant)
allDataMonthly.seas = zeros(height(allDataMonthly),1);
allDataMonthly.seas = month(allDataMonthly.datetime);
allDataMonthly.seas = double(allDataMonthly.seas > 3 & allDataMonthly.seas < 11);

%% Import Spartanburg monthly precip calculate monthly means
[raw] = xlsread('SPGhistPrecip.xlsx');

startDate = datetime(raw(2,1),1,1); %Create monthly datetime vector
endDate = datetime(raw(end,1),12,1);
dataDates = (startDate:calmonths(1):endDate)';

raw = raw(69:end-1,2:end); %Strip away first row and first column to leave only monthly values, use only data from 1950 on
raw(raw == -9999) = NaN; %Replace NoData values with NaN
raw = raw .* 25.4; %Convert to mm

%Calculate monthly average
SPGmonMeanPrecip = nanmean(raw)';
SPGmonSDPrecip = std(raw)';

%% Import Spartanburg monthly precip, calculate monthly mean, convert to PET
[raw] = readtable('SPGhistTemp.csv','ReadVariableNames',0);
raw = raw.Var1; %Convert to just a cell array to facilitate reshaping
raw(end+1:end+5) = {0}; %Pad with extra zeros 
raw = reshape(raw,14,[]);
raw = raw';

raw = raw(68:end-1,:); %Extract only data from 1950 on
raw = cellfun(@str2num,raw(:,2:end-1)); %Convert to numbers and remove first and last columns
raw = (raw - 32).*5/9; %Convert to farenheit

%Calculate monthly average
SPGmonMeanTemp = nanmean(raw);
SPGmonSDTemp = std(raw);

%Convert to PET
load dayLength

daylength = daylength * 24;
monthDays = [30 31 30 31 31 28 31 30 31 30 31 31]';

SPGhistHeatIndex = sum((SPGmonMeanTemp./5).^1.514);
SPGhistAlpha = 6.75e-7*SPGhistHeatIndex^3 - 7.71e-5*SPGhistHeatIndex^2 + 1.792e-2*SPGhistHeatIndex + 0.49239;
SPGhistPET = 16.*(daylength./12).*(monthDays./30).*(10.*SPGmonMeanTemp'./SPGhistHeatIndex).^SPGhistAlpha; %mm/month
SPGhistPETsd = 16.*(daylength./12).*(monthDays./30).*(10.*SPGmonSDTemp'./SPGhistHeatIndex).^SPGhistAlpha; %mm/month

%% Calculate annual wetness index (P/PET)
WY2015WI = sum(allDataMonthly.precip(WY2015I))/sum(allDataMonthly.PET(WY2015I));
WY2016WI = sum(allDataMonthly.precip(WY2016I))/sum(allDataMonthly.PET(WY2016I));
WY2017WI = sum(allDataMonthly.precip(WY2017I))/sum(allDataMonthly.PET(WY2017I));

%% Calculate annual RR 
WY2015RR = sum(allDataDaily.runoff(WY2015I))/sum(allDataDaily.precip(WY2015I));
WY2016RR = sum(allDataDaily.runoff(WY2016I))/sum(allDataDaily.precip(WY2016I));
WY2017RR = sum(allDataDaily.runoff(WY2017I))/sum(allDataDaily.precip(WY2017I));

%% Exceedance analysis w/ hourly/daily data (Flow Duration Curves)

%Hourly
% totalFDC = [sort(allDataHourly.runoff,'descend') linspace(0,1,length(allDataHourly.runoff))']; %Create FDC matrix (height of data, 2): [Sorted runoffs, 0:1]
% 
% WY2015FDC = [sort(allDataHourly.runoff(WY2015I),'descend') linspace(0,1,length(allDataHourly.runoff(WY2015I)))']; %Same as above but for each water year
% WY2016FDC = [sort(allDataHourly.runoff(WY2016I),'descend') linspace(0,1,length(allDataHourly.runoff(WY2016I)))'];
% WY2017FDC = [sort(allDataHourly.runoff(WY2017I),'descend') linspace(0,1,length(allDataHourly.runoff(WY2017I)))'];

%Daily
totalFDC = [sort(allDataDaily.runoff,'descend') linspace(0,1,length(allDataDaily.runoff))']; %Create FDC matrix (height of data, 2): [Sorted runoffs, 0:1]

WY2015FDC = [sort(allDataDaily.runoff(WY2015I),'descend') linspace(0,1,length(allDataDaily.runoff(WY2015I)))']; %Same as above but for each water year
WY2016FDC = [sort(allDataDaily.runoff(WY2016I),'descend') linspace(0,1,length(allDataDaily.runoff(WY2016I)))'];
WY2017FDC = [sort(allDataDaily.runoff(WY2017I),'descend') linspace(0,1,length(allDataDaily.runoff(WY2017I)))'];

%Percent of flow within 5% (quick flow)
%WY2015
qFlow15i = WY2015FDC(:,2) <= 0.3;
qFlow15 = sum(WY2015FDC(qFlow15i,1))/sum(WY2015FDC(:,1));
%WY2016
qFlow16i = WY2016FDC(:,2) <= 0.3;
qFlow16 = sum(WY2016FDC(qFlow16i,1))/sum(WY2016FDC(:,1));
%WY2017
qFlow17i = WY2017FDC(:,2) <= 0.3;
qFlow17 = sum(WY2017FDC(qFlow17i,1))/sum(WY2017FDC(:,1));
%WYtotal
qFlowtoti = totalFDC(:,2) <= 0.3;
qFlowtot = sum(totalFDC(qFlowtoti,1))/sum(totalFDC(:,1));

%% Calculate Total and Annual Richards-Baker 
%%RBI = sum(|dQ/dT|)/sum(Q); Summations over same time period

totalRBI = sum(abs(diff(allDataHourly.runoff)))/sum(allDataHourly.runoff);

WY2015RBI = sum(abs(diff(allDataHourly.runoff(WY2015I))))...
    /sum(allDataHourly.runoff(WY2015I));
WY2016RBI = sum(abs(diff(allDataHourly.runoff(WY2016I))))...
    /sum(allDataHourly.runoff(WY2016I));
WY2017RBI = sum(abs(diff(allDataHourly.runoff(WY2017I))))...
    /sum(allDataHourly.runoff(WY2017I));

%% Calculate Monthly RBI
diffRHourly = timetable(allDataHourly.datetime,[0; abs(diff(allDataHourly.runoff))]); % |dQ/dT| for RBI from above
pathRMonthly = retime(diffRHourly,'monthly','sum'); %Calculates pathlength, numerator of RBI from above, for each month
pathRMonthly.Properties.VariableNames = {'RBIPath'};
allDataMonthly.RBIPath = pathRMonthly.RBIPath; 
allDataMonthly.RBI = allDataMonthly.RBIPath./allDataMonthly.runoff;
allDataMonthly.RBI(isnan(allDataMonthly.RBI)) = 0;

%% Save dataset
save PrecipAndQData.mat allDataMonthly allDataDaily allDataHourly %Save monthly, daily, hourly
