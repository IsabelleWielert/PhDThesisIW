function get_raw_data

[fname_v, pname_v] = uigetfile('Choose evaluated speed data', '*.txt');
cd(pname_v);
[fname_r, pname_r] = uigetfile('Choose raw data', '*.txt');

v = load([pname_v,fname_v]);
num_ds = v(end,1);
tr = load([pname_r,fname_r]);

for i=1:num_ds
    ind_ds_v = v(:,1) == i;
    ind_ds = tr(:,1) == i;
    tr_n_v = v(ind_ds_v,2);
    t_tr_v = v(ind_ds_v,3);
    
    tr_max = tr_n_v(end);
    for j=1:tr_max
        ind_tr_v = find(tr_n_v == j,1);
        if ~isempty(ind_tr_v);
            ind_tr_v = tr_n_v == j;
            ind_tr = tr(:,5);
            t_tr_j = t_tr_v(ind_tr_v);
            tr_j = tr(ind_tr, :);
            ind_raw = tr_j(:,4) == t_tr_j;
            tr_out(ind_raw,:);
        end
    end
end


