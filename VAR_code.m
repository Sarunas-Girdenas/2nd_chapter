% VAR model code, used in introduction


%Specify the model

% Add needed toolboxes, link to obtain this toolbox is at the bottom of this script

addpath('Codes\VARToolbox')
addpath('Codes\FigureToolbox')
addpath('Codes\FigureToolbox\ExportFig')

% Set font to make figures (PDFs) look nice

set(0,'DefaultTextFontName','Palatino')
set(0,'DefaultAxesFontName','Palatino')

%Load Data

[Y label] = xlsread('var_data_3.xls','Data')
unit=0;

%Labels

labels={'Borrowing','GDP','Interest Rate','Unemployment','Inflation'};
%Since we know that our data is evenly distributed from 1971 to 2012 
%we can create column and call it 'Years'
Years=1971:1:2012;

%Detrend the data

Output=detrend(Y(:,2)/10000); %Detrended output
b=detrend(Y(:,1)); % Detrend borrowing
u=detrend(Y(:,4)); % Detrend unemployment
infl=Y(:,5);

%Put them again into matrix

Y1=[b Output Y(:,3) u infl];


%Estimate VAR

VARout = VARmodel(Y1,1);
VARprint(VARout,labels);

% Compute IRFs and error bands

[IRF IRF_opt] = VARir(VARout,42,'bq');
[IRFinf IRFsup IRFmed] = VARirband(VARout,IRF_opt);
VARirplot(IRFmed,0,labels,'irf',IRFinf,IRFsup)

%Output from above commands is given in PDF files. We also need it to
%compute IRFS


% We want to plot IRFs of borrowing shock

%Change Years to Periods

Periods=1:1:length(Years);

%Borrowing
figure;
subplot(2,2,1)
plot(Periods,-IRF(:,1,1)) % Plotting IRF itself agains years
hold;
plot(Periods,-IRFsup(:,1,1),'k') %plot the lower bound
plot(Periods,-IRFinf(:,1,1),'r') %plot the upper bound
hline=refline(0) %add line at zero
title('Borrowing')
xlabel('Periods');
ylabel('Percentage Deviation');

%Output
subplot(2,2,2)
plot(Periods,-IRF(:,2,1)) % Plotting IRF itself agains years
hold;
plot(Periods,-IRFsup(:,2,1),'k') %plot the lower bound
plot(Periods,-IRFinf(:,2,1),'r') %plot the upper bound
hline=refline(0) %add line at zero
title('Output')
xlabel('Periods');
ylabel('Percentage Deviation');

%Interest Rate
subplot(2,2,3)
plot(Periods,-IRF(:,3,1)) % Plotting IRF itself agains years
hold;
plot(Periods,-IRFsup(:,3,1),'k') %plot the lower bound
plot(Periods,-IRFinf(:,3,1),'r') %plot the upper bound
hline=refline(0) %add line at zero
title('Interest Rate')
xlabel('Periods');
ylabel('Percentage Deviation');

%Unemployment
subplot(2,2,4)
plot(Periods,-IRF(:,4,1)) % Plotting IRF itself agains years
hold;
plot(Periods,-IRFsup(:,4,1),'k') %plot the lower bound
plot(Periods,-IRFinf(:,4,1),'r') %plot the upper bound
hline=refline(0) %add line at zero
title('Unemployment')
xlabel('Periods');
ylabel('Percentage Deviation');

%IRFs to GDP equation
%Borrowing
figure;
subplot(2,2,1)
plot(Periods,IRF(:,1,2)) % Plotting IRF itself agains years
hold;
plot(Periods,IRFsup(:,1,2),'k') %plot the lower bound
plot(Periods,IRFinf(:,1,2),'r') %plot the upper bound
hline=refline(0) %add line at zero
title('Borrowing')
xlabel('Periods');
ylabel('Percentage Deviation');

%Output
subplot(2,2,2)
plot(Periods,IRF(:,2,2)) % Plotting IRF itself agains years
hold;
plot(Periods,IRFsup(:,2,2),'k') %plot the lower bound
plot(Periods,IRFinf(:,2,2),'r') %plot the upper bound
hline=refline(0) %add line at zero
title('Output')
xlabel('Periods');
ylabel('Percentage Deviation');

%Interest Rate
subplot(2,2,3)
plot(Periods,IRF(:,3,2)) % Plotting IRF itself agains years
hold;
plot(Periods,IRFsup(:,3,2),'k') %plot the lower bound
plot(Periods,IRFinf(:,3,2),'r') %plot the upper bound
hline=refline(0) %add line at zero
title('Interest Rate')
xlabel('Periods');
ylabel('Percentage Deviation');

%Unemployment
subplot(2,2,4)
plot(Periods,IRF(:,4,2)) % Plotting IRF itself agains years
hold;
plot(Periods,IRFsup(:,4,2),'k') %plot the lower bound
plot(Periods,IRFinf(:,4,2),'r') %plot the upper bound
hline=refline(0) %add line at zero
title('Unemployment')
xlabel('Periods');
ylabel('Percentage Deviation');


%Some remarks about matrix ordering. Notice, that IRF matrix has 3
%dimensions. (a,b,c) where a denotes IRFs, b denotes variable, and c
%denotes shock. So IRFs of first variable to 1st shock is given by (:,1,1).

%And we are done!
%Files needed to run this code could be found at https://sites.google.com/site/ambropo/examples