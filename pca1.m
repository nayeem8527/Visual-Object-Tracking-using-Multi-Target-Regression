function [P] = pca1(feats,d)
[P,~,~] = PCA(feats,d);
% X = feats;
% meanmat = repmat(mean(X,2),1,size(X,2));
% X_norm=X-meanmat;
% [m, ~] = size(X_norm);
% sigma = (X_norm'*X_norm)/m;
% [U,S,V] = svd(sigma);
% P = X_norm*U(:,1:180);

% P=P';
P = normc(P);
