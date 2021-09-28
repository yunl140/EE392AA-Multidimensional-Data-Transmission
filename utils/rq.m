function [R,Q] = rq(A, t)
%rq Triangular-orthogonal decomposition
%   [R,Q] = rq(A), where A is m-by-n, produces an m-by-n upper triangular
%   matrix R and an n-by-n unitary matrix Q so that A = R*Q'.    
%   
%   [R,Q] = rq(A,0) produces the "economy size" decomposigion.
%   If m<n, only the last m columns of R and Q are computed. If m>=n, this
%   is teh same as [R,Q] = rq(A).

reverse_idx_func = @(X) X(end:-1:1, end:-1:1);

if nargin == 1
    [std_Q, std_R] = qr(reverse_idx_func(A)');
else
    [std_Q, std_R] = qr(reverse_idx_func(A)', 0);
end
    Q = reverse_idx_func(std_Q);
    R = reverse_idx_func(std_R)';
end

