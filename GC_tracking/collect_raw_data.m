function collect_raw_data

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

tr_all = zeros(1000000,5);
num_tr = 0;

for f=1:num_folders
    folder = [char(folder_path(f)),'\'];
    folder_results = dir([folder, 'results*']);
    file_name = folder_results.name(9:end);
    file_tr_mat = ['track_', file_name, '.mat'];

    data = load([folder, folder_results.name, '\', file_tr_mat]);            
    tr = data.tr;      
    data = load([folder, folder_results.name, '\data_all.txt']); 
    tr_id = data(:,1);
    
  %data = load([folder, folder_results.name, '\Data_Speed_all_versus_time.mat']); 
  %tr_id = data.data_all(:,3);
    
    
    l_tr = tr_id(end);
    for i=1:l_tr
        id = find(tr_id == i, 1);
        if ~isempty(id)
            num_tr = num_tr + 1;
            ind_tr = find(tr(:,4) == i);
            null_id = find(tr_all(:,1)==0,1);
            ind = null_id:(null_id + length(ind_tr) - 1);
            tr_all(ind, 1) = f;
            tr_all(ind, 2:5) = tr(ind_tr, :);
        end            
    end
    waitbar(f/num_folders);
end
close(h)
tr_all(find(tr_all(:,1)==0,1):end,:)=[];
[FileName,PathName] = uiputfile('*.txt');

path_data = [PathName, FileName]; 
path_num_tr = [PathName, FileName(1:end-4), '_num_tracks.txt']; 
save(path_data , 'tr_all', '-ascii', '-double' ,'-tabs'); 
save(path_num_tr , 'num_tr', '-ascii'); 

assignin('base', 'tr_all', tr_all)