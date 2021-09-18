function [R0,R1,R2] = calc_R0(Net, alpha, mu, gamma, lambda, kappa, N)
adj = Net{1,5}{1,1};
% all(NetRnd{1,5}{1,1} == NetRnd{1,5}{1,2})
% mu(1);       % S_s  -> S_a
% mu(2);       % S_a  -> S_s
% gamma(1);    % I1_a -> S_a
% gamma(2);    % I2_a -> S_a
% lambda(1);   % I1_s -> I1_a
% lambda(2);   % I1_a -> I1_s
% kappa(1);    % I2_s -> I2_a
% kappa(2);    % I2_a -> I2_s

% % old formula as it was in original capstone
% I = eye(N);
% muuuu = mu(1)/(mu(1) + mu(2));
% R1 = (alpha(1)*adj*muuuu + lambda(1)*I)/(gamma(1)+lambda(2));
% R2 = (alpha(2)*adj*muuuu + kappa(1)*I)/(gamma(2)+kappa(2));

% new derivation as of June 2021
R1 = (alpha(1)*adj*mu(1))/(mu(1)+mu(2))/(gamma(1)+lambda(2));
R2 = (alpha(2)*adj*mu(1))/(mu(1)+mu(2))/(gamma(2)+kappa(2));

R1 = EIG2(R1);
R2 = EIG2(R2);
% R1 = max(eigs(R1));
% R2 = max(eigs(R2));
R0 = max(R1, R2);


