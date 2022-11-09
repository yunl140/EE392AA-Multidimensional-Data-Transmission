H = [1, 0.8; 0.6, 1];
U = size(H,1);
G = H.^2;
bvec = [0.4,0.3]';
w = [1,1]';
err = 1e-5;

bvec = 2*log(2)*bvec;
Eu = zeros(U,1);
theta = zeros(U);
A = zeros(U,U,U);
errvec = ones(U,1);
for u = 1:U
    [A(:,:,u), theta(:,u)] = startEllipse(reshape(H(u,:), 1, U, 1), bvec, w);
end
while 1
    for u = 1:U % receiver U
        % Rxx step
        [stheta, order] = sort(theta(:,u), 'descend'); % everything after u in order
        dec_after = find(order==u);
        Itmp = 1 + G(u,order(1:dec_after-1))*Eu(order(1:dec_after-1),1);
        for idx = dec_after:U
            Eu(order(idx)) = (exp(bvec(order(idx))) - 1) * Itmp / G(u,order(idx));
            Itmp = Itmp + Eu(order(idx))*G(u,order(idx));
        end
        % theta step
        for v = 1:U % calculate updated rates at receiver v
            if u == v
                continue;
            end
            [stheta, order] = sort(theta(:,v), 'descend');
            dec_after = find(order==v);
            Stmp = G(v, order)'.*Eu(order);
            Itmp = cumsum([1;Stmp]);
            btmp = log(Itmp(2:end)./Itmp(1:end-1));
            bu = bvec;
            bu(order(dec_after:end)) = btmp(dec_after:end);
            g = bu - bvec;
            errvec(v) = sqrt(g'*A(:,:,v)*g);
        end
    end
    break;
end