function warping_path_w_DGD(t, x)
    
    function printMessage(hObject, eventdata, button, curr)
                    
                    if strcmp(button, 'yes')
                        fprintf('%f,1\n', curr);
                        
                    else
                        fprintf('%f,0\n', curr);
                    end
                    close(gcf);
                end
    
    % Define a list of list
    
    % Iterate over the main list
    for i = 1:length(x)
        % Create a figure with a plot
        x0 = t(x{i}{2}:x{i}{2} - 1 +14);
        x1 = t(x{i}{3}:x{i}{3} - 1 +14);
        
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
            


        % Add buttons to the figure
        uicontrol('Style', 'pushbutton', 'String', 'Yes', 'Position', [20 20 50 20], 'Callback', {@printMessage, 'yes', x{i}{4}});
        uicontrol('Style', 'pushbutton', 'String', 'No', 'Position', [90 20 50 20], 'Callback', {@printMessage, 'no', x{i}{4}});
        
        % Wait for the user to press a button
        waitfor(gcf);
    end
end