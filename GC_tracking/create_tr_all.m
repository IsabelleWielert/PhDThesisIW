function tr_all = create_tr_all(mode)

path = uigetdir(cd, 'Choose folder containing track data');
cd(path)

files = dir('*.mat');
fnames = {files.name}';
n_f = size(files,1);
tr_all = zeros(1000000, 5);
n_ds = 0;
for i=1:n_f
    data = load(fnames{i});
    if strcmp(mode, 'high')
        if isfield(data, 'tr_high')
            tr = data.tr_high;
        else
            continue;
        end
    elseif strcmp(mode, 'low')
        if isfield(data, 'tr_low')
            tr = data.tr_low;
        else
            continue;
        end
    end
    rows = size(tr,1);
    null_id = find(tr_all(:,1)==0,1);
    ind = null_id:null_id + rows - 1;
    tr_all(ind,1) = tr(:,1) + n_ds;
    tr_all(ind, 2:end) = tr(:,2:end);
    n_ds = n_ds + tr(end,1);
end

tr_all(find(tr_all(:,1)==0,1):end,:)=[];
