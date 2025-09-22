function index_set = ARP(V)
[n,r] = size(V);

Vk = V; index_set = [];
for k = 1:r
    % New leverage scores
    pj = vecnorm(Vk(:, k:r), 2, 2).^2 / (r - k + 1);

    % Sample
    new_ind = randsample(n, 1, 'true', pj);
    index_set = [index_set new_ind];

    % Householder matrix
    Qk = zeros(r, r);
    if k > 1
        Qk(1:(k-1), 1:(k-1)) = eye(k-1);
    end

    % Make HH reflector
    x = Vk(new_ind, k:r)';
    e1 = zeros(size(x)); e1(1) = 1;
    u = x - norm(x)*e1; 
    u = u / norm(u);

    % Finish HH matrix
    Qk(k:end, k:end) = eye(r - k + 1) - 2 * (u * u');

    % Eliminate
    Vk = Vk * Qk;
end

end