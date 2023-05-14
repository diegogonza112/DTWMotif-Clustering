function k_driver(symptom)
r=readtable('VMQ_data.xls');
for in = 105:143
    part = r(r{:, "PartID"} == in, :);
    sympt = part{:,  symptom};
    clean = sympt(~isnan(sympt(:,1)),:);
    helper(clean, 14, 4, in);
end
end

function [out] = find_cluster(v, o, run)
x = 1;
    while v(4, x) < 1.35
        o{end + 1} = [v(2, x), v(3, x), v(4, x), run];
        x = x + 1;
    end
    if length(o) > 1
    out = o;
    
    else
    out={};
    end
end

function helper(ts, l, w, z)

[~,~,~, val] = DTWMotifDiscoveryGUI(ts,l,w);
top_k1 = val(2, 1);
[out] = find_cluster(val, [], 0);
ts(top_k1:top_k1+l - 1) = 0 + (120).*rand(l,1);
top_k2 = val(1, 1);
ts(top_k2:top_k2+l - 1) = 0 + (120).*rand(l,1);

run_num = 1;
while val(4, 1) < 1.35
    [~,~,~, val] = DTWMotifDiscoveryGUI(ts,l,w);
    top_k1 = val(2, 1);
    [out] = find_cluster(val, out, run_num);
    ts(top_k1:top_k1+l - 1) = 0 + (120).*rand(l,1);
    run_num = run_num + 1;
end

x = length(out);
if x >= 1
    fil = fopen("BipolarMood.txt", "a");
    for i=1:x
        fprintf(fil, '{%d, %d, %f} Cluster %d from Participant %d\n', out{i}(1), out{i}(2), out{i}(3), out{i}(4), z);
    end
    fclose(fil);
else
    fprintf('None Found for %d\n', z);
end
end