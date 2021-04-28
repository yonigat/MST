%%
% Sanity check for Hydrogen
% Read Hydrogen XS txt  file
e_vec = 0.1:0.1:300;
N = 6.022e23;
Hmat = csvread('Hydrogen.csv');
Hmat = Hmat(:,1:end-1);

%%
res = zeros(length(e_vec),1);
for i = 1:length(e_vec)
    res(i) = sum_sigma_calc(e_vec(i));
end
%%

loglog(e_vec, res)
hold on
loglog(e_vec, Hmat(2,:))
legend({'computed','simulation'})

% figure()
% loglog(e_vec, Hmat(2,:)./res')


% Hmat(2,i)
% sum_sigma
% (Hmat(2,i) - sum_sigma)/Hmat(2,i)


