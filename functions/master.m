% from input and num_fr, writes master file
% input is xls tracking file from fiji
function [M,num_tr] = master(input)

    global num_fr DT trl_t;
%       i = 1; j = 1; M = NaN((max(input(:,1))+1)-min(input(:,1)),num_fr,6); % x,y,vx,vy,w,quality
       i = 1; j = 1; M = NaN((max(input(:,2))+1),num_fr,6); % x,y,vx,vy,w,quality
%       i = 1; j = 1; M = NaN(40000,num_fr,4); % x,y,vx,vy,w,quality
    while 1
        foo=find(input(:,2)==j-1);
        if length(input(foo,4)) > trl_t
             M(i,input(foo,8)+1,[1:4 6]) = [input(foo,4) input(foo,5) [NaN;diff(input(foo,4))/DT] [NaN;diff(input(foo,5))/DT] input(foo,3)];
   %              M(i,input(foo,8)+1,[1:4 4]) = [input(foo,4) input(foo,5) [NaN;diff(input(foo,4))/DT] [NaN;diff(input(foo,5))/DT]];% positions,velocities, and quality
            i = i+1;
        end
        j = j+1;
        if j == max(input(:,2))+1
            break;
        end
    end
    num_tr = i-1; clear i j% new number of tracks

    M = M(1:num_tr,:,:); % trim MASTER to new num_tr
%     figure(); scatter(M(:,10,1),M(:,10,2),'.'); 
end

 