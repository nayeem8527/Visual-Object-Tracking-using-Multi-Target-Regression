function [pred,GT] = classification(data,dict,gt)

W = dict.D;
S = dict.W;
% X = data;
X = data';
n1sq = sum(X.^2,1);
n1 = size(X,2);
D = (ones(n1,1)*n1sq)' + ones(n1,1)*n1sq -2*X'*X;
K = exp(-D/(2*0.01^2));
X = K*X';

% W = A*X';
Z = W*X;
pred = S*Z;
GT = S*Z;
Y = gt;
Y = repmat(Y,size(pred));
Y = Y(1:2,:);
% pred = sqrt(sum((pred-Y).^2,1)/2);

% sparse coding
% Dict = dict.D; % 200 X 496
% 
% thres=0.001;
% s=Dict*data;% added this line
% st=@(s,threshold)sign(s).*max(0,abs(s)-threshold);%added this line
% Gamma=st(s,thres);% added this line
% 
% % reconstruction error
% TW = dict.TW / sum(dict.TW);
% avgT = dict.T * TW';
% d = pinv(dict.D) * Gamma;
% A_norm = sqrt(sum(d .* d));
% DX = d ./ (ones(size(d,1),1) * A_norm + eps);
% Y = repmat(avgT, 1, size(data,2));
% dif1 = Y - DX;
% rec_error = sqrt(sum(dif1.^2));
% % classification error
% H = repmat([1;0],1,size(data,2));
% dif2 = H - dict.W * Gamma;
% cls_error = sqrt(sum(dif2.^2));

pred = 0.8*pred(1,:) + 0.2*pred(2,:);
end