% function [Rxx, bsum , bsum_lin] = SWF(Eu, h, user_ind, Rnn, N , cb)
%
% the inputs are: 
% Eu  the U x 1 energy vector. If a single scalar, same energy for all users
%     Each user energy should be scaled by N/(N+nu)if there is cyclic prefix 
%     This energy is the trace of the corresponding user Rxx (u)
% h   The TIME-DOMAIN Ly x sum(Lx(u)) x N channel for all users
% user_ind  The start index for each user, in the same order as Eu
% Rnn The Ly x Ly noise autocorrelation
% N   The number of used tones (equally spaced over (0,1/T) at N/T.
% cb  cb = 1 for complex, cb=2 for real baseband
%
% the outputs are:
% Rxx A block-diagonal psd matrix with the input autocorrelation for each
%     user on each tone. Rxx has size (sum(Lx(u)) x sum(Lx(u)) x N .
%     sum trace(Rxx) over tones and spatial dimensions equal the Eu 
% bsum the maximum rate sum.
% bsum  bsum_lin - the maximum sum rate with a linear receiver
% 
%     b is an internal convergence (vector, rms) value, but not sum rate

function [Rxx, bsum, bsum_lin] = SWF(Eu, h, user_ind, Rnn, N , cb)
U = numel(user_ind);
if numel(Eu) == 1
    Eu = Eu*ones(U);
end

i = 0;
[Ly, Lx_sum, taps] = size(h);
Rxx = diag(zeros(Lx_sum, N));
b = zeros(U, 1);

Lxu=zeros(1,U);
user_ind=[ user_ind Lx_sum + 1];
for u=1:U
    Lxu(u)=user_ind(u+1) - user_ind(u);
end


H = 1/sqrt(N)*fft(h, N, 3);

for u = 1:U
   
    st = user_ind(u);
    
    if u < U
        en = user_ind(u+1)-1;
    else
        en = Lx_sum;
    end
    
    for n = 1:N
        Rxx(st:en, st:en,n) = N*Eu(u)*eye(en-st+1);
    end
end

b = zeros(U);
while 1
    b_prev = b;
    
    for u = 1:U
        %new user
        st = user_ind(u);

            if u < U
                en = user_ind(u+1)-1;
            else
                en = Lx_sum;
            end
        M_u = zeros(en-st+1,en-st+1, N);
        g_u = [];
        for n = 1:N
            %new tone
            Rnn_un = Rnn;
            for v = 1:U
                if v ~= u
                    st = user_ind(v);
                    if v < U
                        en = user_ind(v+1)-1;
                    else
                        en = Lx_sum;
                    end
                    Hvn = H(:, st:en, n);
                    Rxxvn = Rxx(st:en, st:en, n);
                    Rnn_un = Rnn_un + Hvn*Rxxvn*Hvn';
                end
            end
            
            st = user_ind(u);

            if u < U
                en = user_ind(u+1)-1;
            else
                en = Lx_sum;
            end
            
            
            Hun = H(:, st:en, n);
            
            Hun_til = sqrtm(inv(Rnn_un))*Hun;

            [F, D, M_n] = svd(Hun_til);
            s = svd(Hun_til);
            
            M_u(:,:,n) = M_n;
            g_u(:,n) = s.^2;
           
        end
        
        g_flat = reshape(g_u, [numel(g_u),1]);
        [g_sort, ind] = sort(g_flat, 'descend');
        
        
        L = numel(g_sort);
        E_bar = N*Eu(u);
        j = L;
        
        while 1
            
            K_til = E_bar+sum(1./g_sort(1:j));
            K = K_til/j;
            
            
            Ej = K - 1/g_sort(j);
            
            if Ej <0
                j = j-1;
            else               
                E = K-1./g_sort(1:j);
                B=.5*log2(K*g_sort(1:j));
                L_star = j;
                break
            end
            
        end
        
        e = zeros(L,1);
        e(ind(1:L_star))=E;
        
        e = reshape(e, [en-st+1, N]);
        
   
        b(u,1) = sum(B);
        
        for n = 1:N
            Rxx(st:en, st:en, n) = M_u(:,:,n)*diag(e(:,n))*M_u(:,:,n)';
        end
    end

    i = i+1;
    
    if norm(rms(b_prev-b)) <= 1e-5
        break
    end
    if i>1000
        break
    end
end
%b = b(:,1);
%Hcell = mat2cell(H,Ly,Lx_sum,ones(1,N));
%Hexpand = [blkdiag(Hcell{1,1,:})]
%Rcell = mat2cell(Rxx,Ly,Lx_sum,ones(1,N));
%Rcell = mat2cell(Rxx,Lx_sum,Lx_sum,ones(1,N));
%Rxxexpand = [blkdiag(Rcell{1,1,:})]
%bsum = log2(det(eye(Lx_sum*N)+Hexpand*Rxxexpand*Hexpand'));
%bsum = log2(det(eye(Ly*N)+Hexpand*Rxxexpand*Hexpand'));
bsum=0;
for n=1:N
    bsum=bsum+(1/cb)*log2(det(eye(Ly)+H(:,:,n)*Rxx(:,:,n)*H(:,:,n)'));
end
bsum=real(bsum);
bs=zeros(1,U);
bsum_lin=0;
for u=1:U
    indices=user_ind(u):user_ind(u)+Lxu(u)-1;
    for n=1:N
     bs(u)=bs(u)+(1/cb)*(log2(det(eye(Ly)+H(:,:,n)*Rxx(:, ...
       :,n)*H(:,:,n)')) - log2(det(eye(Ly)+H(:,:,n)*Rxx(:, ...
       :,n)*H(:,:,n)'- H(:,indices,n)*Rxx(indices, ...
       indices,n)*H(:,indices,n)')));
    end
    bsum_lin=bsum_lin+real(bs(u));
end
Rxx=Rxx/N;
