function [E, theta, b] = minPMAC(H, bu, w)
% This is the ellipsoid method part and the main function to call.

% the inputs are: 
% 1)  H, an Ly by U by N channel matrix. Ly is the number of receiver antennas, 
%     U is the total number of users and N is the total number of tones.
%     H(:,:,n) is the channel matrix for all users on tone n
%     and H(:,u,n) is the channel for user u on tone n. In this code we assume each user 
%     only has single transmit antenna, thus H(:,u,n) is a column vector. 
% 2)  w, a U by 1 vector containing the weights for each user's power.
% 3)  bu, a U by 1 vector containing the target rates for all the users.

% the outputs are:
% 1)  E, a U by N matrix containing the powers for all users and over all tones
%     that minimizes the weighted-sum power. E(u,:) is the power allocation for
%     user u over all tones.
% 2)  theta, the optimal U by 1 dual variable vector containing optimal weights
%     of rates. theta determines the decoding order.
% 3)  b, a U by N matrix containing the rates of all users on all tones after convergence.


err = 1e-9;                              % error tolerance
count = 0;                               
[Ly, U, N] = size(H);

[A, g] = startEllipse(H, bu, w);         % starting ellipsoid 
theta = g;

bu = bu  * log(2);                       % conversion from bits to nuts 

while 1
                                         % Ellipsoid method starts here
    [f, b, E] = Lag_dual_f(H, theta, w, bu);
    g = sum(b,2) - bu;                   % sub-gradient
    
    if sqrt(g' * A * g) <= err           % stopping criteria
        break
    end
                                         % Updating the ellipsoid
    tmp = A*g / sqrt(g' * A * g);
    theta = theta - 1 / (U + 1) * tmp;
    
    A = U^2 / (U^2 - 1) * (A - 2 / (U + 1) * (tmp * tmp'));
    
    ind = find(theta < zeros(U,1));
    
    while ~isempty(ind)                  % This part is to make sure that theta is feasible,
        g = zeros(U,1);                  % it was not covered in the lecture notes and you may skip this part
        g(ind(1)) = -1;
        tmp = A * g / sqrt(g' * A * g);
        theta = theta - 1 / (U + 1) * tmp;
        A = U^2 / (U^2 - 1) * (A - 2 / (U + 1) * (tmp * tmp'));
        ind = find(theta < zeros(U,1));
    end   
    count = count+1;
end

b = b /log(2);                            % conversion from nuts to bits
   
    
    
        

    