function eval_data_tr

%global switching window in units of minutes
min_t = 19.5;
max_t = 21;


%open dialog for choosing txt.file
[FileName,PathName] = uigetfile('*.txt',...
                           'Choose text file containing posistion data!');
%if cancel is pushed, then return.
if (isscalar(FileName) == 1) && (isscalar(PathName) == 1);
    return;
end  
cd(PathName)
tr_all = load([PathName, FileName]);
num_ds = tr_all(end,1);


%dissect data into high mode, low mode and global switching window
rows = size(tr_all,1);
tr_gs = zeros(rows,5);
tr_high = zeros(rows,5);
tr_low = zeros(rows,5);
num_tr_gs = 0;
num_tr_high = 0;
num_tr_low = 0;

for i=1:num_ds
    ind_s = tr_all(:,1)==i;
    tr_s = tr_all(ind_s,2:end);
    last_tr = tr_s(end,1);
    h = waitbar(0,['Dataset: ', int2str(i),'. Dissect data ...']);
    for j=1:last_tr
        tr_id = find(tr_s(:,1)==j);
        if ~isempty(tr_id)
            %data inside global switching window
            time_tr = tr_s(tr_id, 2);
            ind_t = find(time_tr>min_t & time_tr<max_t);
            if ~isempty(ind_t)
                num_tr_gs = num_tr_gs + 1;
                dp = length(ind_t);
                null_ind = find(tr_gs==0,1);
                ind = null_ind:null_ind+dp-1;
                tr_gs(ind,1) = i;
                tr_gs(ind,2:end) = tr_s(ind_t + tr_id(1) - 1, :);
            end
            %data of high velocity mode
            ind_t = find(time_tr < min_t);
            if ~isempty(ind_t)
                num_tr_high = num_tr_high + 1;
                dp = length(ind_t);
                null_ind = find(tr_high(:,1)==0,1);
                ind = null_ind:null_ind+dp-1;
                tr_high(ind,1) = i;
                tr_high(ind,2:end) = tr_s(ind_t + tr_id(1) - 1,:);
            end
            %data of low velocity mode
            ind_t = find(time_tr > max_t);
            if ~isempty(ind_t)
                num_tr_low = num_tr_low + 1;
                dp = length(ind_t);
                null_ind = find(tr_low(:,1)==0,1);
                ind = null_ind:null_ind+dp-1;
                tr_low(ind,1) = i;
                tr_low(ind,2:end) = tr_s(ind_t + tr_id(1) - 1,:);
            end
            waitbar(j/last_tr);
        end
    end
    close(h)
end

tr_gs(find(tr_gs(:,1)==0,1):end,:)=[];
tr_high(find(tr_high(:,1)==0,1):end,:)=[];
tr_low(find(tr_low(:,1)==0,1):end,:)=[];

if ~isempty(tr_gs)
    assignin('base', 'tr_gs', tr_gs);
    assignin('base', 'num_tr_gs', num_tr_gs);
    path_gs = [PathName, 'tr_gs.txt'];
    save(path_gs , 'tr_gs', '-ascii', '-double' ,'-tabs'); 
    path_gs = [PathName, 'num_tr_gs.txt'];
    save(path_gs , 'num_tr_gs', '-ascii', '-double' ,'-tabs'); 
end
if ~isempty(tr_high)
    assignin('base', 'tr_high', tr_high);
    assignin('base', 'num_tr_high', num_tr_high);
    path_high = [PathName, 'tr_high.txt'];
    save(path_high , 'tr_high', '-ascii', '-double' ,'-tabs'); 
    path_high = [PathName, 'num_tr_high.txt'];
    save(path_high , 'num_tr_high', '-ascii', '-double' ,'-tabs'); 
end
if ~isempty(tr_low)
    assignin('base', 'tr_low', tr_low);
    assignin('base', 'num_tr_low', num_tr_low);
    path_low = [PathName, 'tr_low.txt'];
    save(path_low , 'tr_low', '-ascii', '-double' ,'-tabs'); 
    path_low = [PathName, 'num_tr_low.txt'];
    save(path_low , 'num_tr_low', '-ascii', '-double' ,'-tabs'); 
end