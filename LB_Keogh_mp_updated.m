function [mp,mpi] = LB_Keogh_mp_updated(ts, subseqlen, minlag, warpmax, DoNotCompute)
% Computes all pairs LB_Keogh in a 1 to many fashion
% Query is used to generate the envelope in each case
% Inputs are time series: ts,
% subsequence length (integer): subseqlen
% minimum index offset between subsequence pairs (integer): minlag
% maximum warping offset (integer):  warpmax
% DoNotCompute: vector of true and false booleans showing which subsequences to
% prune

% Outputs are matrix profile and index using LB_Keogh

if ~isvector(ts) || ~isvector(DoNotCompute)
    error('expected a vector input');
end

istransposed = isrow(ts);
if istransposed
    ts = transpose(ts);
end
if isrow(DoNotCompute)
    DoNotCompute = transpose(DoNotCompute);
end
if subseqlen < 4 || minlag < 4
    error('bad parameters');
end

subcount = length(ts) - subseqlen + 1;
mu = movmean(ts, [0 subseqlen-1], 'Endpoints', 'discard');
sig = movstd(ts, [0 subseqlen-1],  1, 'Endpoints', 'discard');

subs = zeros(subseqlen, subcount);


% Normalize everything in advance. If this is memory prohibitive,
% subdivide your problem and do the same thing. Just don't normalize every
% single time on the fly. That is generally slower and more complicated.
for i = 1 : subcount
    subs(:, i) = (ts(i : i + subseqlen - 1) - mu(i)) ./ sig(i);
end

% These are stored in
mp = inf(subcount, 1);
mpi = inf(subcount, 1);

U = NaN(subseqlen, subcount);
L = NaN(subseqlen, subcount);


if warpmax > 0
    for i = 1 : subcount
        if DoNotCompute(i)
            continue;
        end
        U(:, i) = movmax(subs(:, i), [warpmax warpmax]);
        L(:, i) = movmin(subs(:, i), [warpmax warpmax]);
    end
else
    for i = 1 : subcount
        if DoNotCompute(i)
            continue;
        end
        U(:, i) = subs(:, i);
        L(:, i) = subs(:, i);
    end
end

% In C++ code, you would set this up in a slightly different manner. It
% always takes some memory to set all of these up, but you can store 1
% sided U and L in O(n) memory, which can be used to quickly form
% per subsequence U,L with simple merges. Next, you can tile this.
% You may repeat the normalization step a few times, but each time you
% normalize and form U,L for some subset of the problem, you use them to form
% O(tilelen^2) pairs. This way the cost of data movement is trivialized at
% typical subsequence lengths relative to the cost of the main loop of each tile. 


for i = 1 : subcount
    for j = i + minlag : subcount        
        if DoNotCompute(i) || DoNotCompute(j)
            continue;
        end
        d = max(sum((U(:, i) - subs(:, j)).^2 .* (U(:,i) < subs(:, j)) + (L(:, i) - subs(:, j)).^2 .* (L(:, i) > subs(:, j))),...
                sum((U(:, j) - subs(:, i)).^2 .* (U(:,j) < subs(:, i)) + (L(:, j) - subs(:, i)).^2 .* (L(:, j) > subs(:, i))));       
        if d < mp(i)
            mp(i) = d;
            mpi(i) = j;
        end
        if d < mp(j)
            mp(j) = d;
            mpi(j) = i;
        end
    end
end

mp = sqrt(max(0, mp));

if istransposed
    mp = transpose(mp);
    mpi = transpose(mpi);
else
    
end


end

