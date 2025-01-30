function num_tr = find_num_tracks_in_t_int(data, t_l, t_h)

ds = data(:,1);
tr = data(:,2);
t = data(:,3);
ds_max = ds(end);
num_tr = 0;

for i=1:ds_max
    ind_ds = find(ds==i);
    if ~isempty(ind_ds)
        tr_ds = tr(ind_ds);
        tr_max = tr(ind_ds(end));
        for j=1:tr_max
            ind_tr = find(tr_ds==j);
            if ~isempty(ind_tr)
                t_ds = t(ind_ds);
                t_tr = t_ds(ind_tr);
                ind_t = find(t_tr>=t_l & t_tr<= t_h, 1);
                if ~isempty(ind_t)
                    num_tr = num_tr + 1;
                end
            end
        end
    end
end
