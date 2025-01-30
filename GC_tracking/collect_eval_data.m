function collect_eval_data

%Choose Folder containing movie folders
Main_Path = uigetdir(cd, 'Choose folder containing movie folders');
cd(Main_Path);
folder_path = cell(100,1);
Main_Path_dir = dir(Main_Path);
k = 1;
for z=1:size(Main_Path_dir,1)
   if Main_Path_dir(z).isdir == 1 && strcmp(Main_Path_dir(z).name, '.')==0 ...
       && strcmp(Main_Path_dir(z).name, '..')==0; 
      folder_path(k) = cellstr([Main_Path, '\', Main_Path_dir(z).name]);
      k = k + 1;
   end
end
folder_path(k:end)=[];
num_folders = size(folder_path, 1);

h = waitbar(0,'Please wait...');
data_all = zeros(1000000,6);

num_tr = 0;
for f=1:num_folders
    folder = [char(folder_path(f)),'\'];
    folder_results = dir([folder, 'results*']);
    %cd([folder, folder_results.name, '\']);
    %file_name = folder_results.name(9:end);
    %file_tr_mat = ['track_', file_name, '.mat'];
    %file_tr_txt = ['track_', file_name, '.txt'];

    %data = load([folder, folder_results.name, '\', file_tr_mat]);            
    %data = data.tr;   
    %tr = load([folder, folder_results.name, '\', file_tr_txt]);

    %data_all = load([folder, folder_results.name, '\data_all_Speed_versus_time.txt']); 
    %tr_id = data_all(:,3);
    data = load([folder, folder_results.name, '\data_all.txt']); 
    tr_id = data(:,1);
    length(tr_id)
    
    l_tr_id = tr_id(end);
    num_dead_tr = 0;
    for i=1:l_tr_id
        id = find(tr_id == i, 1);
        if isempty(id)
            num_dead_tr = num_dead_tr + 1;
        end
    end
    null_id = find(data_all(:,1)==0,1);
    ind = null_id:(null_id + length(tr_id) - 1);
    data_all(ind, 1) = f;
    data_all(ind, 2:6) = data;
    num_tr = num_tr + (l_tr_id - num_dead_tr);
    waitbar(f/num_folders);
end
close(h)
data_all(find(data_all(:,1)==0,1):end,:)=[];

[FileName,PathName] = uiputfile('*.txt');
drawnow
path_data = [PathName, FileName]; 
path_num_tr = [PathName, FileName(1:end-4), '_num_tracks.txt']; 
save(path_data , 'data_all', '-ascii', '-double' ,'-tabs'); 
save(path_num_tr , 'num_tr', '-ascii'); 

assignin('base', 'data_all', data_all)