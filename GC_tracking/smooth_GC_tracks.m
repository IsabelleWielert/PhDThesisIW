function smooth_GC_tracks(data, tr_nice)

win_sm = 50;

tr = data(:,1);
t = data(:,2);
v = data(:,5);
i_max = length(tr_nice);
for i=1:i_max
    ind2 = tr==tr_nice(i);
    t_name = ['t_', int2str(tr_nice(i))];
    t_tr = t(ind2);
    assignin('base', t_name, t_tr)
    v_name = ['v_sm_', int2str(tr_nice(i))];
    v_sm_tr = smooth(v(ind2), win_sm);
    assignin('base', v_name, v_sm_tr)
    plot(t_tr, v_sm_tr, '-k')
    hold on
end
