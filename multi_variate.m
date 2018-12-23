clear
clc

load X
load Y

Y_new = zeros(2,size(Y,2));
Y_new(1,1:500)=1;
Y_new(2,501:size(Y,2))=1;
Y=Y_new;

% lambda = 0.00001; 
% beta = 0.001;
% gamma = 0.01;

eta = 0.00001;
K = X'*X;
K_inv = pinv(K);
N = size(X,2);

% vec = -5:1:3;
% val = ones(9,3);
% temp=vec';
% val(:,1) = val(:,1).*temp;
% val(:,2) = val(:,2).*temp;
% val(:,3) = val(:,3).*temp;
% combs = nchoosek(vec,3);
% combs = [combs;val];
% results = zeros(size(combs,1),10);
% for j=1:size(combs,1)
%     disp(j);
lambda = 10^(-5); 
beta = 10^(-4); 
gamma = 10^(-3); 
S = rand(size(Y,1),size(Y,1));
for i=1:10

    temp1 = S'*S;
    temp2 = (lambda.*(N.*K_inv));
    temp3 = ((S'*Y)*K_inv);
    A = sylvester(temp1,temp2,temp3);

    [U,SS,V] = svd(S);
    SS_inv = pinv(SS);
    G_S = ((-2/N)*((Y-((S*A)*K))*(A*K)')) + ((((beta.*U)*SS_inv)*abs(SS))*V') + ((2*gamma).*S);
    S = S - (eta.*G_S);

    % for checking the convergence
    error = (1/N)*trace((Y-((S*A)*K))'*(Y-((S*A)*K))) + lambda*trace((A*K)*A') + beta*(trace(sqrt(S'*S))) + gamma*(trace(S'*S));
    disp(error);    
end
% end

% W = A*X';
% Z = W*X;
% Y_pred = S*Z;

% -5,-4,-3 or -4,-4,-4