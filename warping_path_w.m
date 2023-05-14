function [cost, path] = warping_path_w(s1, s2, w)

x0 = s1;
x1 = s2;

n = length(x0);
m = length(x1);

if exist('w', 'var')
    w = max(w, abs(m-n)); % adapt window size (*)
    DTW = inf(n,m);
    for i = 2:n
        for j = max(2, i-w):min(m, i+w)
            DTW(i, j) = 0;
        end
    end
    DTW(1,1) = 0;

    for i = 2:n
        j_min = max(2, i-w);
        j_max = min(m, i+w);
        for j = j_min:j_max
            cost_arr = abs(x0(i) - x1(j));
            DTW(i, j) = cost_arr + min( [DTW(i-1, j) DTW(i, j-1) DTW(i-1, j-1)] );
        end
    end

else

    DTW = zeros(n,m);
    DTW(1,:) = inf;
    DTW(:,1) = inf;
    DTW(1,1) = 0;

    for i0 = 2:length(x0)
        for i1 = 2:length(x1)
            cost_arr = abs(x0(i0) - x1(i1));
            DTW(i0, i1) = cost_arr + min( [DTW(i0-1, i1) DTW(i0, i1-1) DTW(i0-1, i1-1)] );
        end
    end

end

[cost_arr, path] = min(DTW, [], 2);
cost = cost_arr(end);

end