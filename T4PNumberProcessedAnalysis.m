    %% This program collect all Data from T4P_Contour, determines Pilus production rate and the nummber of present pili
    clear all;
    close all;
    clc;

    %% read in files:
    datafiles{1,1}= dir(fullfile('path of your folder with .mat data from T4P_Contour analysis'));
    framerate=19.333; % Adjust when you used another framerate
    Frame_ana= 150; % decided to only analyse the first 7.5s due to phototoxicity


    %% Determine Present pili and pilus production rate
    PresentPilidistri=[];
    Npilipminute=[];


    for z=1:length(datafiles)

        datafilesi=datafiles{1,z}



        for ii = 1:length(datafilesi)
            tmp = load(fullfile(datafilesi(ii).folder, datafilesi(ii).name));

            datastr = sprintf('Contour_%u', ii);  % Generate data string
            data.(datastr) = tmp;
        end




        B=cell2mat(struct2cell(data));



        PresentPili=[];
        NumberPili=[];
        Npilipframe=[];
        Npilipminute=[]
        for i=1:length(B)


            Mat=B(i).Pili;


            Mat=B(i).Pili;
            [L,n] = bwlabel(Mat',8);
            Lneu=L;
            % Close line if there are less than 10 frames between two
            % detected pili, pre-processing matrix
            SE = strel('line',10,0)
            Lneu=imclose(Lneu, SE);
            [L,n] = bwlabel(Lneu,8);
            Pii=zeros(size(L));
            for b=1:n
                [row,col]=find(L==b);
                rowcol=[row col];
                rowstar=rowcol(1,1);
                if size(rowcol,1)>10
                    for y=rowcol(:,2)
                        Pii(rowstar,y)=1;
                    end
                end
            end
            Mat=Pii;
            Mat=imclose(Pii, SE);
            [L,n]=bwlabel(Mat(1:end,:));

            F=imbinarize(L);
            Mat=F';

            % count present pili for each frame:
            Pili_present=[];
            for g=1:Frame_ana
                Pilithere= sum(Mat(g,:));
                Pili_present=[Pili_present; Pilithere];
            end
            PresentPili=[PresentPili; Pili_present];

            % Loop to make sure, that pili are connected in the right way
            % and that there are at least 10 frames break if pili should be
            % counted, separately, and loop for counting if a new pilus
            % appears in matrix

            einsen=1;
            newpilipertime=[];
            for p=4:size(Mat,2)-7
                countp=0;
                for k = 2:Frame_ana
                    if k>6
                        if Mat(k-5,p)==0 && Mat(k,p)==0 && Mat(k+1,p)==einsen && Mat(k+4,p)==einsen && Mat(k+6,p)==einsen &&Mat(k+5,p)==einsen &&Mat(k+7,p)==einsen &&Mat(k+8,p)==einsen &&Mat(k+9,p)==einsen
                            countp=countp+1;
                        end
                    else
                        if Mat(k+1,p)==einsen && Mat(k,p)==0 && Mat(k+4,p)==einsen && Mat(k+6,p)==einsen &&Mat(k+5,p)==einsen
                            countp=countp+1;
                        end
                    end
                end
                newpilipertime=[newpilipertime;countp];
            end

            % Calculate the new number of pili appearing for each frame
            NEWpiliperframe=sum(newpilipertime)/Frame_ana;
            NEWpiliperminute=framerate*60*NEWpiliperframe; % Number of pili per minute

            Npilipframe=[Npilipframe;NEWpiliperframe];
            Npilipminute=[Npilipminute;NEWpiliperminute];


        end


        meannumber=fitdist(PresentPili,'Normal');
        PresentPilidistri{z}=PresentPili;
        NumberParams(z,1:2)=meannumber.Params;

        meanpiliperminute=fitdist(Npilipminute,'Normal');
        Npilipminutedistri{z}=Npilipminute;
        PiliperminuteParam(z,1:2)=meanpiliperminute.Params;
    end
    %%  Analysis of Pilus number

    figure();
    h=histogram(PresentPilidistri{1,1},'Normalization','probability')
    h.FaceColor=[0.8 0.8 0.8]
    h.EdgeColor=[0.8 0.8 0.8]
    hold on
    ylabel('Probability')
    xlabel('Present pili')
    title(' Mean number of present pili')
    % filenamefigure = [[savepath,'\'],'Pilihisto.png'];
    % saveas(gcf, filenamefigure);
    hold off


    figure();
    x=[PresentPilidistri{1,1}]
    g1 = repmat({'Data'},length(PresentPilidistri{1,1}),1);
    g = [g1];
    boxplot(x,g, 'symbol', '')
    yt = get(gca, 'YTick');
    axis([xlim    floor(min(yt))-0.5  ceil(max(yt))+5.5])
    xt = get(gca, 'XTick');
    hold on
    ylabel('Number of pili')
    title(' Mean number of present pili')
    xticklabels('Data')
    xtickangle(45)
    hold off
    % you want to save the figure?
    % filenamefigure = [[savepath,'\'],'Presentpili.png'];
    % saveas(gcf, filenamefigure);

    %% Analysis of Productionrate [T4P/Minute]

    figure();
    x=[Npilipminutedistri{1,1}]
    g1 = repmat({'Data'},length( Npilipminutedistri{1,1}),1);
    g = [g1];
    boxplot(x,g, 'symbol', '')
    yt = get(gca, 'YTick');
    axis([xlim    floor(min(yt))-0.5  ceil(max(yt))+5.5])
    xt = get(gca, 'XTick');
    hold on
    ylabel('Pili per minute')
    title('Number of pili per minute')
    xticklabels('Data')
    xtickangle(45)
    hold off
    %     filenamefigure = [[savepath,'\'],'Npilidistri.png']; % If you want to
    %     save the figure, uncomment
    %     saveas(gcf, filenamefigure);
    figure();

    h=histogram(Npilipminutedistri{1,1},'Normalization','probability')
    h.FaceColor=[0.8 0.8 0.8]
    h.EdgeColor=[0.8 0.8 0.8]
    hold on
    legend('data')
    ylabel('Probability')
    xlabel('Pili per minute')
    title('Pili per minute')
    % filenamefigure = [[savepath,'\'],'Npilihisto.png'];
    % saveas(gcf, filenamefigure);
    hold off




    %% Statistical tests if you want to analyse more than one conditon and compare

    % ksdauer=[];
    % for i = 1:2
    % [l p]=kstest2(Dauerdistri{1,3},Dauerdistri{1,i});
    % ksdauer=[ksdauer;l p];
    % end
    %
    % kselo=[];
    % for i = 1:2
    % [l p]=kstest2(Elongationdistri{1,3},Elongationdistri{1,i});
    % kselo=[kselo; l p];
    % end

