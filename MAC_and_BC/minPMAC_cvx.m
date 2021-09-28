function [Eun, theta, bun] = minPMAC_cvx(H, bu_min, w)
% minPMAC_cvx solves minimize sum_u(w_u*tr(Rxx(u))), subject to
% bu>=bu_min
% This function uses CVX package

% the inputs are: 
% 1)  H, an Ly by U by N channel matrix. Ly is the number of receiver antennas, 
%     U is the total number of users and N is the total number of tones.
%     H(:,:,n) is the channel matrix for all users on tone n
%     and H(:,u,n) is the channel for user u on tone n. In this code we assume each user 
%     only has single transmit antenna, thus H(:,u,n) is a column vector. 
% 2)  bu, a U by 1 vector containing the target rates for all the users.
% 3)  w, a U by 1 vector containing the weights for each user's power.

% the outputs are:
% 1)  Eu, a U by N matrix containing the powers for all users and over all tones
%     that minimizes the weighted-sum power. E(u,:) is the power allocation for
%     user u over all tones.
% 2)  theta, the optimal U by 1 dual variable vector containing optimal weights
%     of rates. theta determines the decoding order.
% 3)  bu, a U by N matrix containing the rates of all users on all tones after convergence.


[Ly, U, N] = size(H);
w = reshape(w,U,1);
bu_min = reshape(bu_min,U,1);

err_threshold = 1e-6;
D = eye(U)+diag(-ones(U-1,1),1);
% initialize
err = 1;
theta = 0.01*ones(1,U);
[theta, order] = sort(theta, 'descend');
step_size_init = 1/min(bu_min);
count = 1;
old_f = -1;
bun = zeros(U,N);
Eun = zeros(U,N);
while (err > err_threshold)
    % Rxx step
    for n = 1:N
        cvx_begin quiet
            variable En(U) nonnegative
            expression Ru(U)
            for v = 1:U
                Ru(v) = log_det(H(:,order(1:v))*diag(En(order(1:v)))*H(:,order(1:v))'+eye(Ly));
            end
            Ru = log2(exp(1))*Ru;
            minimize (w'*En - Ru'*pos(D*theta'))
        cvx_end
        for v = 1:U
            Ru(v) = log_det(H(:,order(1:v))*diag(En(order(1:v)))*H(:,order(1:v))'+eye(Ly));
        end
        Ru = log2(exp(1))*Ru;
        bun(order(2:end),n) = diff(Ru);
        bun(order(1),n) = Ru(1);
        Eun(order,n) = En;
    end
    % theta step
    bu = sum(bun,2)
    theta = max(theta + step_size_init/count*(bu_min - bu)', 0);
    [theta, order] = sort(theta, 'descend');
    new_f = w'*En;
    err = abs(new_f - old_f);
    old_f = new_f;
    count = count + 1;
end

[~,rev_order] = sort(order);
theta = theta(rev_order);

end
    
  