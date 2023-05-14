function [pr_rate,best_so_far,motiffirst,motifsec,dtw_time,top_k] = dtw_mpGUI_tweaked(ts, mp_ind, mp_sorted, subseqlen, minlag, warpmax, dnc, best_so_far,gui_name,motiffirst,motifsec,top_k_num)

mu = movmean(ts, [0 subseqlen-1], 'Endpoints', 'discard');
sig = movstd(ts, [0 subseqlen-1], 1, 'Endpoints', 'discard');
subcount = length(ts) - subseqlen + 1;

lb_kim_time = 0;
lb_keogh_time = 0;
dtw_time = 0;
lb_kim_iter = 0;
lb_keogh_iter = 0;
dtw_iter = 0;
ea_cnt = 0;

% top_k = zeros(1,3);
% top_k_num = 20;         %Select the k top motifs

tr= (dnc==1);
%[~,idx] =intersect(mp_ind,tr,'stable');
mp_ind(tr) = -1;
%mp_sorted(idx) = -1;
tic
for i = 1 : subcount      
    id = (mp_ind(i));    
    if id==-1 || dnc(id)==1
        continue;
    end  
    neighbors = [id + minlag : subcount];  
    mp_ind(dnc==1) = -1;
    neigh = mp_ind(id);
    neighbors(neighbors==neigh)= [];
    neighbors = [neigh,neighbors];
    
    for j = neighbors       
        idp = j;        
        if idp==-1 || dnc(idp)==1
            continue;
        end       
        
        a = (ts(id : id + subseqlen - 1) - mu(id))./sig(id);
        b = (ts(idp : idp + subseqlen - 1) - mu(idp))./sig(idp);
        
        %%%%%%%%%%%%%% Compute LB_Kim and skip DTW computation if higher than best-so-far %%%%%%%%%%%%%%%%%%  
        tic
        lb_Kim = max(abs(a(1)-b(1)),abs(a(end)-b(end)));
        temp = toc;        
        lb_kim_time = temp + lb_kim_time;
        if lb_Kim >= best_so_far     
            lb_kim_iter = lb_kim_iter + 1;
             continue;
             
        %%%%%%%%%%%%%% Compute LB_Keogh and skip DTW computation if higher than best-so-far %%%%%%%%%%%%%%%%%%  
        else
            Ua = movmax(a,[warpmax warpmax]);
            La = movmin(a,[warpmax warpmax]);
            tic           
            LB_Keogh = lb_upd(b,Ua,La,best_so_far);            
            temp = toc;
            lb_keogh_time = temp + lb_keogh_time;            
            if sqrt(LB_Keogh) >=best_so_far
                lb_keogh_iter = lb_keogh_iter + 1;
                continue;
            end
        end        
        
        %%%%%%%%%%%%%% Compute the DTW distance if the previous two checks failed %%%%%%%%%%%%%%%%%% 
        tic
        
         %415 731|115 191|215 232|119 58|466 103|65 126

        [dist,ea] = dtw_upd(a,b,warpmax,best_so_far);

%       [dist2,ix,iy] = dtw(a,b,warpmax,'euclidean'); %get a warping path using Signal Toolbox dtw(); use the L2 metric

        % Since we are neglecting best_so_far from here on in,...
        % ...plug in to dtw() routine with ability to return warp_path, and
        % don't bother with early_abandonment (ea) value
        [dist3,warp_path] = warping_path_w(a,b,warpmax);

        temp = toc;
%         ea_cnt = ea + ea_cnt;
        dtw_iter = dtw_iter + 1;
        dtw_time = temp + dtw_time;
        
        %%%%%%%%%%%%%% Update best-so-far if needed %%%%%%%%%%%%%%%%%%
        if dist < best_so_far
%             fprintf('%f \n', best_so_far);
          
            if (~exist('top_k','var') || isempty(top_k))
                top_k = [id idp dist]; 
            else
                top_k = [top_k; [id idp dist]]; 
            end

%           best_so_far = dist;
            motiffirst = id;
            motifsec = idp;     
            gui_name.plotData(id, idp)
%             gui_name.plotMotifdtw(id, idp)
            gui_name.plotMotifdtw(id, idp,warp_path)    %DK tweak
            gui_name.drawgui()
%             pause(1);                         
%             dnc(mp_ind(mp_sorted>best_so_far)>0) = 1;            
%             mp_ind(mp_sorted>best_so_far) = -1;
%             mp_sorted(mp_sorted>best_so_far) = -1;         
        end        
    end
end

%Extract unique elements of top_k[]
top_k_ordered = sortrows([sort(top_k(:,1:2),2) top_k(:,3)],3);  %Order motif indices in each row from lowest to highest
top_k_unique = unique(top_k_ordered,'rows');                    %Then, eliminate duplicate rows
top_k_unique = sortrows(top_k_unique,3);                        %Then, sort by MP, lowest first
top_k_num = min(top_k_num,size(top_k_unique,1));                %If there's less than top_k_num, adjust top_k_num
top_k = top_k_unique(1:top_k_num,:);                            %Finally, take the top_k_num

phase2_time = toc;
% lb_kim_time
% lb_kim_iter
% lb_keogh_time
% lb_keogh_iter
% dtw_time
% phase2_time
% dtw_iter
% ea_cnt
pr_rate = dtw_iter;
end

