function [num_tr] = calc_num_tr(data)

num_tr = 0;

ds_max = data(end,1);

for i=1:ds_max
    ds_id = data(:,1)==i;
    if isempty(find(ds_id, 1))
        continue;
    end
    tr = data(ds_id, 2);
    tr_max = tr(end);
    for j=1:tr_max 
        tr_id = find(tr == j, 1);
        if ~isempty(tr_id)
            num_tr = num_tr +1;
        end
    end
end