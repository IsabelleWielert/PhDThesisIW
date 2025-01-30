    
function plot_speed_hist(speed_all, num_tracks, num_dead, datapoints)
    figure(100),
    
    [n,xout] = hist(speed_all, 0:1/10:5);
    subplot(2,1,1),
    bar(xout,n, .8, 'r'), 
    %fitting a gaussian to speed histogram
    f = fittype('gauss1');
    res_fit  = fit(xout',n',f);
    confi_speed = confint(res_fit);
    err_speed = confi_speed(2,2) - confi_speed(1,2);
    title(['Speed: Histogram (\Delta t = ',...
        num2str(1.0, '%2.1fs)')] , 'FontSize', 14),
    hold on
    plot(res_fit, '-k'),
    hold off
    xlabel('Speed [µm/s]', 'FontSize', 12), ...
    ylabel('Frequency', 'FontSize', 12),
    xlim([0 5]);
    grid on;
    text(2.5,max(n)-50,...
    ['Gaussian Modell:  Y = a*exp(-((x-b)/c)^2).\newline\Rightarrow Av. speed=b=(',...
        num2str(res_fit.b1,'%1.3f'),'\pm', num2str(err_speed, ...
        '%1.3f)µm/s'),'\newline# of tracks: ',num2str(num_tracks, '%2.0f'), ...
        '\newline# of dead tracks: ',num2str(num_dead, '%2.0f'),' (', ...
        num2str(num_dead/num_tracks*100, '%2.0f'),...
        '%)\newline# of datapoints: ', ...
        num2str(sum(datapoints),'%3.0f')], 'FontSize',12),
    subplot(2,1,2),
    [n2,xout2] = hist(datapoints, 1:30:600);
    bar(xout2,n2, .8, 'r'),
    title('Datapoints per track: Histogram' , 'FontSize', 14),
    xlabel('# of datapoints per track', 'FontSize', 12), ...
    ylabel('Frequency', 'FontSize', 12),
    grid on;
end