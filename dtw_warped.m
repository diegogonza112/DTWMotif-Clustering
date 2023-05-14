function dtw_warped(t, s1, s2)
    
    x0 = t(s1:s1+14);
    x1 = t(s2:s2+14);
    
    warp=1:1:size(x0);
    
    figure
    plot(warp, x0, warp, x1);
    hold on
    DTW = zeros(length(x0), length(x1));
    DTW(1,:) = inf;
    DTW(:,1) = inf;
    DTW(1,1) = 0;
    
    for i0 = 2:length(x0)
        for i1 = 2:length(x1)
            cost = abs(x0(i0) - x1(i1));
            DTW(i0, i1) = cost + min( [DTW(i0-1, i1) DTW(i0, i1-1) DTW(i0-1, i1-1)] );
        end
    end
    
    [cost, path] = min(DTW, [], 2);
    plot(warp, x1(path));
    legend({'x_0', 'x_1', 'x_1 warped to x_0'});
end