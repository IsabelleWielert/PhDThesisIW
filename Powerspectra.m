%% Program to determine the stiffness via the Powerspectrum and calibrate the trap


close all;
clc;
clear all;
% read in the data

allpath='Your folder with measurements for a certain laser intensity';
% arrays to save the stiffness
K1=[];
K2=[];
K3=[];
K4=[];

% read out data
liste=dir(allpath);

Paths={liste.name};
Paths=Paths(:,3:end);
B={};
for o=1:length(Paths)
    
    pfadnamen=fullfile(allpath,Paths{o})
    B{o}=pfadnamen;
end



for o=1:length(B)
    path=B{o};
    
    d = path;
    
    files = dir([d '/*.txt']);
    % Absolute Dateinamen erstellen
    files = strcat([d '/'], {files.name});
    
    for s=1:2
        for h=1:2
            data = textread(files{s});
            data=data';
            
            
            Masstab=0.886525*10^(-7);
            framerate=2500;
            
            
            %Average of data
            
            c=0;
            L_data=2500*5;
            max_length=25000;
            FFT_mittel=zeros(1250*5+1,1);
            % Fourier transformation of data to get power spectral density
            for i=1000:L_data:max_length+1000
                
                Fx=fft(data(h,i:i+L_data)-mean(data(h,i:i+L_data)));
                Fx = Fx(1:floor(end/2)+1);
                
                f=0:1/(1/framerate)/(length(Fx)-1):framerate;
                
                PSD = (1/(framerate*length(Fx))).*(Fx.*conj(Fx));
                
                PSD(2:end-1) = 2*PSD(2:end-1);
                
                
                
                for j=1:6251
                    FFT_mittel(j)= FFT_mittel(j)+PSD(j);
                    
                end
                c=c+1;
            end
            
            PSD=FFT_mittel./c;
            
            PSD=PSD';
            %fit intervall
            
            
            drag=6*pi*0.5*10^(-6)*10^(-3)*0.891;
            
            % in x-direction
            
            [xData, yData] = prepareCurveData(f, PSD);
            
            % Set up fittype and options. To get cut off frequency with which the
            % stiffness can be determined, beacuse the spectral density should follow:
            %Pk=D/(2*pi^2(fc^2+fk^2)) and the stiffness can then be detrmined from the
            %ctu-off frequency fc
            ft = fittype( '(a./(b^2.+x.^2))', 'independent', 'x', 'dependent', 'y' );
            
            excludedPoints = excludedata( xData,yData, 'domain',  [45 300]);
            opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
            opts.Display = 'Off';
            opts.StartPoint = [10,200];
            opts.Exclude = excludedPoints;
            opts.Lower = [1,10];
            opts.Upper = [100,400];
            % Fit model to data5.
            [fitresult, gof] = fit( xData, yData, ft, opts );
            
            % Plot fit with data. Uncomment if you would like to see it
            
            % f1=figure( 'Name', ' Data 1' );
            % figure(f1);
            % z = plot( fitresult, xData, yData, excludedPoints);
            % %h = plot( fitresult, xData, yData);
            %
            % legend( z, 'Data', 'Excluded data', 'Fit', 'Location', 'NorthEast' );
            % %legend( h, 'Data', 'Fit', 'Location', 'NorthEast' );
            % % Label axes
            % xlabel Hz
            % ylabel Powerdensity
            % grid on
            % %
            
            
            coefs=coeffvalues(fitresult);
            
            k_gem=coefs(2)*2*3.141592654*drag; %determines the stiffness via Powerspectrum via k=2pi*friction coefficient*cut-off frequency
            
            % Saves the stiffness
            if h==1 && s==1 %obenx
                K1=[K1 k_gem];
            elseif h==1 && s==2%untenx
                K2=[K2 k_gem];
            elseif h==2 && s==1%obeny
                K3=[K3 k_gem];
            elseif h==2 && s==2%unteny
                K4=[K4 k_gem];
            end
            
            % Figure to make sure that fitting worked:
            
            f2=figure('Position',[300 ,300, 1000, 500],'name',['Powerspectrum gemittelt Data1 mit excluded Peak ' num2str(path)]);
            figure(f2);
            loglog(xData,yData,'.');
            hold on;
            loglog(f,coefs(1)./(coefs(2).^2.+f.^2));
            
            xlabel('Hz');
            ylabel('Powerdensity');
            
            
            legend('data','fit');
        end
    end
end
format short
% average over all data to get mean stiffness for a specific intensity
kx_oben=mean(K1) % kx oben
kx_unten=mean(K2) % kx unten
ky_oben=mean(K3) % ky oben
ky_unten=mean(K4) % ky unten
