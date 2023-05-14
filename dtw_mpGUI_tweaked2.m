function [pr_rate,best_so_far,motiffirst,motifsec,dtw_time,top_k,top_k_num] = dtw_mpGUI_tweaked2(ts, mp_ind, mp_sorted, subseqlen, minlag, warpmax, dnc, best_so_far,gui_name,motiffirst,motifsec,top_k_num)

mu = movmean(ts, [0 subseqlen-1], 'Endpoints', 'discard');
sig = movstd(ts, [0 subseqlen-1], 1, 'Endpoints', 'discard');
subcount = length(ts) - subseqlen + 1;

top_k = zeros(subcount,3); % DK tweak
top_k_idx = 0;              % DK_tweak
mp_ind_tmp = mp_ind;        %DK tweak - stash original copy of mp_ind for later

lb_kim_time = 0;
lb_keogh_time = 0;
dtw_time = 0;
lb_kim_iter = 0;
lb_keogh_iter = 0;
dtw_iter = 0;
ea_cnt = 0;

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
        temp = toc;
        ea_cnt = ea + ea_cnt;
        dtw_iter = dtw_iter + 1;
        dtw_time = temp + dtw_time;

        %%%%%%%%%%%%%% Update best-so-far if needed %%%%%%%%%%%%%%%%%%
        if dist < best_so_far

            best_so_far = dist;
            motiffirst = id;
            motifsec = idp;
            gui_name.plotData(id, idp)
            %             gui_name.plotMotifdtw(id, idp)
            gui_name.plotMotifdtw(id, idp,[])    %DK tweak
            gui_name.drawgui()
            pause(1);
            dnc(mp_ind(mp_sorted>best_so_far)>0) = 1;
            mp_ind(mp_sorted>best_so_far) = -1;
            mp_sorted(mp_sorted>best_so_far) = -1;
        end
    end
end
% phase2_time = toc;

%%%%   DK Kludge Feb 22 '23 - re-compute all DTW distances for the top
%%%%   motif paired with all other points and extract the top-k. Ignore any
%%%%   prior pruning except for Euclidean Distance (which was implemented
%%%%   before dtw_mpGUI_tweaked2() was called.

% Restore mp_ind to its original state (with only ED constraints
mp_ind = mp_ind_tmp;
idp = motiffirst;

for i = 1 : subcount
    id = (mp_ind(i));
    % if id==-1  || dnc(id)==1
    if id==Inf              % DK_tweak
        continue;
    end

    a = (ts(id : id + subseqlen - 1) - mu(id))./sig(id);
    b = (ts(idp : idp + subseqlen - 1) - mu(idp))./sig(idp);

%   [dist,ea] = dtw_upd(a,b,warpmax,best_so_far);
    [dist,~] = dtw_upd(a,b,warpmax,Inf);   % DK tweak - set best_so_far to Inf so dtw_upd won't early-abandon
%   [dist3,warp_path] = warping_path_w(a,b,warpmax); % DK_tweak - get warping path but discard distance

    top_k_idx = top_k_idx + 1;
    top_k(top_k_idx,:) = [id idp dist];

    gui_name.plotData(id, idp)
    gui_name.plotMotifdtw(id, idp, [])          % DK tweak - add null arg
%   gui_name.plotMotifdtw(id, idp,warp_path)    % DK_tweak - plot using warp path
    gui_name.drawgui()


end

pr_rate = dtw_iter;

%Extract unique elements of top_k[]
% top_k_ordered = sortrows([sort(top_k(:,1:2),2) top_k(:,3)],3);  %Order motif indices in each row from lowest to highest
% top_k_unique = unique(top_k_ordered,'rows');                    %Then, eliminate duplicate rows
top_k_unique = unique(top_k,'rows');                            %Eliminate duplicate rows
top_k_unique = sortrows(top_k_unique,3);                        %Then, sort by MP, lowest first

% The following should select the top_k_num members of top_k_unique with lowest distances that
% aren't within warpmax of any current member of top_k

i = 1;  % Skip the null
j = 0;
top_k = zeros(top_k_num,3);

while (i < size(top_k_unique,1))   % Don't take more than top_k_num motifs and stop if all top_k_unique have already been scanned
    i = i + 1;  % Look at the next line in top_k_unique
    k = top_k_unique(i,1); % What index are we matching?

    if (abs(motiffirst-k) > warpmax)    % Don't take anything within warpmax of top motif, motiffirst
        if (all(top_k(:) == 0))          % If there's no content in top_k, accept the first entry
            j = j + 1;
            top_k(j,:) = top_k_unique(i,:);
        elseif (all(abs(k-top_k(1:j,1)) > warpmax)) % Otherwise, don't take anything within warpmax of any prior motif
            j = j + 1;
            top_k(j,:) = top_k_unique(i,:);
        end
    end
end

top_k = top_k(1:j,:);
top_k_num = min(top_k_num,j);        %If there's less than top_k_num, adjust top_k_num

end