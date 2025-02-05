close all; 
clc; 
clear all; 

cmap={[160 160 160]./256; [30 30 30]./256; []./256;[20 149 0]./256};

set(0, 'defaultFigurePosition', [100 100 330 300]);
cmap=[[128 128 128]./256; [77 77 77]./256; [254 128 0]./256;[255 0 0]./256];

set(0,'DefaultAxesColorOrder',cmap)
set(0,'defaultAxesFontSize',10);
set(groot, 'DefaultAxesFontName', 'Arial')
cmap={[128 128 128]./256; [77 77 77]./256; [254 128 0]./256;[255 0 0]./256};
cmap1=[[128 128 128]./256; [77 77 77]./256; [254 128 0]./256;[255 0 0]./256];
%% All raw Data of the countings 

%read in the data for each condition, here it was for the strains: G4, Clp24, Clp17 and
%Clp32 

xG4=[	0	0	0	0	2	2	2	4	4	4	6	6	6	8	8	8	10	10	10	10	10	13	13	13	16	16	16	19	19	19]
yG4=[	40000	42600	93000	74000	11000	33800	34400	480000	442000	102000	180000	940000	490000	1880000	1390000	2000000	2000000	1200000	9700000	1960000	10000000	2640000	3400000	4660000	10000000	16000000	15800000	9600000	43200000	30500000];



x24=[	0	0	0	0	2	2	2	4	4	4	6	6	6	8	8	8	10	10	10	10	10	13	13	13	16	16	16	19	19	19];
y24=[	35000	26600	29400	34400	9200	35000	40800	290000	340000	72000	180000	226000	320000	1650000	1490000	2840000	3740000	2090000	1175000	3460000	4080000	3760000	9400000	4000000	14100000	32400000	29200000	30400000	102000000	156000000];

x17=[0	0	0	0	2	2	2	4	4	4	6	6	6	8	8	8	10	10	10	10	10	13	13	13	16	16	16	19	19	19];
y17=[32000	85000	46000	49000	61000	146000	85000	212000	710000	292000	4340000	4680000	6600000	7400000	21100000	26100000	30000000	32800000	33800000	18750000	43600000	18800000	55000000	48000000	11200000	22000000	31800000	4700000	16000000	23300000];



x32=[0	0	0	0	2	2	2	4	4	4	6	6	6	8	8	8	10	10	10	10	10	13	13	13	16	16	16	19	19	19];
y32=[43000	81000	85000	43000	16300	198000	148000	690000	1450000	472000	4100000	6300000	8500000	20100000	20300000	18800000	35800000	35700000	38000000	50000000	20600000	71000000	43200000	35200000	30000000	42000000	10600000	10600000	23000000	33200000];

xvals={xG4,x24,x17,x32};
yvals={yG4,y24,y17,y32};

%% prepare the fit 
figure();
for i=1:4
    
    logy=log(yvals{i});
    x=xvals{i};
    %% look at raw data
    
    %
    %     scatter(xvals{i},logy)
    %     hold on
    
    %% linear fit in linear regime!
    linx=[]
    liny=[]
    for l =1:length(x)
        
       
        
    for j = [ 2 4 6 8 ] % just to crop the data to linear regime

            
            if j==x(l)
                
                linx=[linx;x(l)];
                liny=[liny;logy(l)];
            end
        end
    end
    
    
    Linearegimex{i}=linx;
    Linearegimey{i}=liny;
    
end


%% statistical tests-follow instructions for ANOVA TEST, Linear Regressions
figure();
groupi=[1 1 1 1 1 1 1 1 1 1 1 1   2 2 2 2 2 2 2 2 2 2 2 2   3 3 3 3 3 3 3 3 3 3 3 3   4 4 4 4 4 4 4 4 4 4 4 4 ]';
xi=cell2mat(Linearegimex)
xii=reshape(xi,[],1)
yi=cell2mat(Linearegimey)
yii=reshape(yi,[],1)
gscatter(xii,yii,groupi,'bgrk','x.o+')




%%
figure(); 
 cars=table(xii,yii,groupi);
  cars.groupi = categorical(cars.groupi)
    cars.xii = categorical(cars.xii)


aoctool([xii(1:13); xii(27:39)],[yii(1:13); yii(27:39)],[groupi(1:13); groupi(27:39)]);
aoctool([xii(1:13); xii(14:26)],[yii(1:13); yii(14:26)],[groupi(1:13); groupi(14:26)]);
aoctool([xii(1:13); xii(40:end)],[yii(1:13); yii(40:end)],[groupi(1:13); groupi(40:end)]);

x1=[xii(1:12); xii(13:24)];
y1=[yii(1:12); yii(13:24)];
g1=[groupi(1:12); groupi(13:24)];

x2=[xii(1:12); xii(25:36)];
y2=[yii(1:12); yii(25:36)];
g2=[groupi(1:12); groupi(25:36)];


x3=[xii(1:12); xii(37:end)];
y3=[yii(1:12); yii(37:end)];
g3=[groupi(1:12); groupi(37:end)];
%% Pooled Data for planktonic and aggegrating strains

pool1=[xii];
pooly1=[yii]; 
poolgroup=[groupi(1:12);groupi(1:12);groupi(13:24);groupi(13:24)];
carsgroup=table(pool1,pooly1,poolgroup);
carsgroup.poolgroup=categorical(carsgroup.poolgroup);
aoctool(pool1,pooly1,poolgroup)
fitpool=fitlm(carsgroup,'pooly1~pool1*poolgroup')
w = [2 4 6 8];
figure()
gscatter(xii,yii,groupi,'kkrr','xxoo')
line(w,feval(fitpool,w,'1'),'Color','k','LineWidth',2)
line(w,feval(fitpool,w,'2'),'Color','r','LineWidth',2)

hpooled=anova(fitpool)
xlabel('Time [h]')
ylabel('log(CFU)')
legend('off')
set(gca,...
    'FontSize',10,...
    'FontWeight','bold',...
    'FontName','Arial')
% set(gca, 'Position', get(gca, 'OuterPosition') - ...
%     get(gca, 'TightInset') * [-1 0 1 0; 0 -1 0 1; 0 0 1 0; 0 0 0 1]);
ax=gca;
ax.LineWidth=1.5
hold off 


%% Pairwise significance analysis, important are the interaction terms * 



cars1=table(x1,y1,g1);
cars2=table(x2,y2,g2);
cars3=table(x3,y3,g3);

 cars.groupi = categorical(cars.groupi);
cars1.g1 = categorical(cars1.g1);
cars2.g2 = categorical(cars2.g2);
cars3.g3 = categorical(cars3.g3);
 
fit = fitlm(cars,'yii~xii*groupi','Intercept','false')

aoctool(xii,yii,groupi);
CM=fit.CoefficientCovariance
SE = diag(sqrt(CM))
fit1 = fitlm(cars1,'y1~x1*g1')
fit2 = fitlm(cars2,'y2~x2*g2')
fit3 = fitlm(cars3,'y3~x3*g3')
% 

CM=fit1.CoefficientCovariance
SE = diag(sqrt(CM))

h=anova(fit)
h1=anova(fit1,'summary')
h2=anova(fit2,'summary')
h3=anova(fit3,'summary')



