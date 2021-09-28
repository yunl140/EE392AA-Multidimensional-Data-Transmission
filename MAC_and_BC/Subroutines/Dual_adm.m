function [f, b, E] = Dual_adm(G, theta, w, bu, Eu);
% this function computes the Lagrange dual function by solving the
% optimization problem (calling the function minPtone) on each tone.
% the inputs are: 
% 1)  G, an Ly by U by N channel matrix. Ly is the number of receiver antennas, 
%     U is the total number of users and N is the total number of tones.
%     G(:,:,n) is the channel matrix for all users on tone n
%     and G(:,u,n) is the channel for user u on tone n. In this code we assume each user 
%     only has single transmit antenna, thus G(:,u,n) is a column vector. 
% 2)  theta, a U by 1 vector containing the weights for the rates.
% 3)  w, a U by 1 vector containing the weights for each user's power.
% 4)  bu, a U by 1 vector containing the target rates for all the users.
% 5)  Eu, a U by 1 vector containing the power constraints for all the
%     users. 

% the outputs are:
% 1)  f, the Lagrange dual function value.
% 2)  b, a U by N vector containing the rates for all users and over all tones
%     that optimizes the Lagrangian. b(u,:) is the rate allocation for user
%     u over all tones.
% 3)  E, a U by N vector containing the powers for all users and over all tones
%     that optimizes the Lagrangian. E(u,:) is the power allocation for
%     user u over all tones.


[Ly, U, N] = size(G);
f = 0; 
b = zeros(U,N);
E = zeros(U,N);

% Performing optimization over all N tones,
for i = 1:N
    [temp, b(:,i), E(:,i)] = minPtone(G(:,:,i), theta, w); 
    f = f + temp;
end

f = f - theta' * bu + w' * Eu;