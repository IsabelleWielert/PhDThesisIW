function [v_all, n_tracks] = create_speed_all(tr_all, fps, delta_t, min_t)

%delta_t in units of seconds. Converted to dist in units of frames
dist = round(fps*delta_t);
%min_t is minimum duration of single track in units of minutes
min_dp = min_t*60*fps;
n_ds = tr_all(end, 1);
v_all = zeros(size(tr_all,1), 4);
n_tracks = 0;

for i=1:n_ds
    ind_ds = find(tr_all(:,1) == i);
    if isempty(ind_ds)
        continue;
    else
        tr_ds = tr_all(ind_ds,2:end);
        n_tr = tr_ds(end,1);
        for j=1:n_tr
            ind_tr = find(tr_ds(:,1)==j);
            if length(ind_tr)>min_dp
                n_tracks = n_tracks + 1;
                t = tr_ds(ind_tr, 2);
                pos = tr_ds(ind_tr, 3:4);
                v = calc_speed(pos, t.*60, dist);
                dp = length(v);
                t_v = t(1:end-dist);
                null_id = find(v_all(:,1)==0,1);
                ind = null_id:(null_id + dp - 1);
                v_all(ind,1) = i;
                v_all(ind,2) = j;
                v_all(ind,3) = t_v;
                v_all(ind,4) = v;
            end
        end
    end
end
v_all(find(v_all==0,1):end,:)=[];

function speed = calc_speed(pos, time, dist)

tr_len = length(time);
dp = tr_len - dist;
speed = zeros(dp,1);
for j=1:dp;
    ind = j + dist;
    s_vec = [pos(ind,1)-pos(j,1);pos(ind,2)-pos(j,2)];
    delta_t = time(ind)-time(j);   
    speed(j) = norm(s_vec)/delta_t;
end