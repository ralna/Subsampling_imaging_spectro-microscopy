function [err2, errF] = compareUpToColPerm(A, B,abs)
    
    if nargin == 2
        abs =0;
    end

    numCol = size(A,2);
    allPerms = perms(1:numCol);

    err2 = Inf; errF = inf;
    for jj = 1:length(allPerms)
        newerr2 = norm(A - B(:, allPerms(jj,:)));
        if newerr2 < err2
            err2 = newerr2;
        end
        newerrF = norm(A - B(:, allPerms(jj,:)), 'fro');
        if newerrF < errF
            errF = newerrF;
        end
    end
    if abs ~=1
    err2 = err2/norm(A,2);
    errF = errF/norm(A,"fro");
    end
end