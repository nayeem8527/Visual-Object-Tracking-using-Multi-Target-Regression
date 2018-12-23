function dict = update_dict(feat, dictpars, dict)

X = feat.feaArr;

% A = dict.D;
S = dict.W;
% Y = feat.boxes;

temp_y = zeros(size(X,2),10);
Y = feat.boxes;
temp_y(:,1:2)= Y(:,1:2);
temp_y(:,3) = Y(:,1)+Y(:,3);
temp_y(:,4) = Y(:,2);
temp_y(:,5) = Y(:,1)+Y(:,3);
temp_y(:,6) = Y(:,2)+Y(:,4);
temp_y(:,7) = Y(:,1);
temp_y(:,8) = Y(:,2)+Y(:,4);
temp_y(:,9) = (Y(:,1)+Y(:,3))/2;
temp_y(:,10) = (Y(:,2)+Y(:,4))/2;
Y = temp_y';

% Y = feat.label;
% newY = zeros(2,size(Y,2));
% ind1 = find(Y==1);
% ind2 = find(Y==-1);
% newY(1,ind1)=1;
% newY(2,ind2)=1;
% Y=newY;
eta = 0.000001;
K = X'*X;
K_inv = pinv(K);
N = size(X,2);

lambda = 10^(-5); 
beta = 10^(-4); 
gamma = 10^(-3); 

for i=1:6000

    temp1 = S'*S;
    temp2 = (lambda.*(N.*K_inv));
    temp3 = ((S'*Y)*K_inv);
    A = sylvester(temp1,temp2,temp3);

    [U,SS,V] = svd(S);
    SS_inv = pinv(SS);
    G_S = ((-2/N)*((Y-((S*A)*K))*(A*K)')) + ((((beta.*U)*SS_inv)*abs(SS))*V') + ((2*gamma).*S);
%     G_S = ((-2/N)*((Y-((S*A)*K))*(A*K)')) + ((((beta.*U)*SS_inv)*abs(SS))*V') + ((gamma).*S);
    
    S = S - (eta.*G_S);

    % for checking the convergence
%     error = (1/N)*trace((Y-((S*A)*K))'*(Y-((S*A)*K))) + lambda*trace((A*K)*A') + beta*(trace(sqrt(S'*S))) + gamma*(trace(S'*S));
%     disp(error);    
end
% error = (1/N)*trace((Y-((S*A)*K))'*(Y-((S*A)*K))) + lambda*trace((A*K)*A') + beta*(trace(sqrt(S'*S))) + gamma*(trace(S'*S));
% disp(error);    
W = A*X';
dict.D = W;
dict.W = S;


% online dictionary learning using stochastic gradient descent algorithm
% Inputs:
%    -para: learning paramters
%       .D: initial dictionary
%       .A: initial transform matrix
%       .W: initial classifier
%       .mu: reguliarization parameter
%       .nu1: reguliarization parameter
%       .nu2: reguliarization parameter
%       .maxIters: iteration number
%       .rho: learning rate parameter
%    -trainingdata: input training des. with size of n X N
%    -HMat: label matrix of training des. with size of m X N
%    -QMat: optimal code matrix of training des. with size of K X N
%
% Outputs:
%    -model: learned paramters
%       .D: learned dictionary
%       .A: learned transform matrix
%       .W: learned classifier
%    -stat:
%       .fobj_avg: average objective function value
%
%  Author: Zhuolin Jiang (zhuolin@umiacs.umd.edu)
%



% trn = feat.feaArr;
% num_bases = dictpars.numBases;
% num_iters = dictpars.maxIters;
% num_images = size(trn, 2);
% numClass = dictpars.numcls; % number of objects
% HMat = zeros(numClass, num_bases);
% HMat(1, feat.label == 1) = 1;
% HMat(2, feat.label == -1) = 1;
% 
% %% compute Q
% dictLabel = [];
% numPerClass = num_bases / numClass;
% for classid = 1:numClass
%     labelvector = zeros(numClass,1);
%     labelvector(classid) = 1;
%     dictLabel = [dictLabel repmat(labelvector,1,numPerClass)];
% end


% QMat = zeros(num_bases, size(trn, 2)); % energy matrix
% for frameid = 1:size(trn, 2)
%     label_training = HMat(:, frameid);
%     [~, maxid1] = max(label_training);
%     for itemid = 1:size(dict.D',2)
%         label_item = dictLabel(:, itemid);
%         [~, maxid2] = max(label_item);
%         if(maxid1 == maxid2)
%             QMat(itemid, frameid) = 1;
%         end
%     end
% end


% gamma = dictpars.gamma; % sparse coding parameters
% lambda = dictpars.lambda;
% nu1 = dictpars.nu1;
% nu2 = dictpars.nu2;
% mu = dictpars.mu;
% rho = dictpars.rho;
% bsize = dictpars.batchSize;
% bsize = num_images / 5;
% n0 = num_images/10
% n0 = num_images/(bsize*10);

% Dpart = dict.D; % dictionary
% W = dict.W; % classifier
% model.A = dict.A; % transform matrix

% param.lambda = dictpars.lambda;
% param.lambda2 = 0;
% param.mode = 2;

% DOING THE SAME WAY AS DONE IN INITIALIZATION
% T0=11;  %sparsity level for each patch
% l2=0.5; %weight for the negative log-determinant penalty in problem formulation
% l4=l2;  %weight for the incoherence penalty in problem formulation
% l3=l2;  %weight for the Frobenius norm penalty used within CG (in transform update step)
% p=1;  %parameter for the incoherence penalty in problem formulation
% mu=1e-9; %step size within CG in the transform update step
% numt=5;  %number of iterations of alternating transform learning algorithm (i.e., iterations of sparse coding and transform update)
% numg=10; %number of iterations of CG in the transform update step
% Data = trn;
% STY = T0*ones(1,size(Data,2));
% [K,n]=size(Dpart); 
% ix=find(STY>0); q=Data(:,ix); STY=STY(:,ix); N=size(q,2);
% ez=K*(0:(N-1));STY=STY + ez;
% ZX=Data*Data';
% for i=1:numt    
%     %%%%%%Sparse Coding Step%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     X1=Dpart*q;
%     [s]=sort(abs(X1),'descend');
%     X = X1.*(bsxfun(@ge,abs(X1),s(STY)));
%     %%%%%%Transform Update Step%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     RSu=(q*X')';
%     %CG iterations
%     for j=1:numg
%         ZZ=(Dpart*Dpart').^(p-1);
%         cc = (-2*RSu) + (2*Dpart*ZX');
%         temp = Dpart'*Dpart;
%         temp = 0.2*eye(size(temp)) + temp;
%         cc2=(-2)*Dpart*(inv(temp));            
%         cc3= 2*Dpart;
%         cc4=2*p*((ZZ*Dpart) - (((diag(ZZ))*(ones(1,n))).*Dpart));
%         deD= l2*cc2 + cc + l3*cc3  + l4*cc4;
%         if(j==1)
%             g=deD;
%             d=-g;
%         else
% %             be = ((norm(deD,'fro'))^2)/((norm(g,'fro'))^2);
%             be = ((kpNorm(deD,min(size(deD)),1))^2)/((kpNorm(g,min(size(g)),1))^2);
%             g=deD;
%             d= - g + be*d;
%         end
%         Dpart = Dpart + (mu)*d;   %constant step size
%     end
%     %post-normalization of the rows of the transform
%     SL=diag(1./sqrt(sum((Dpart').^2)));
%     Dpart=SL*Dpart;
% end
% % param1.D = Dinit;
% Dinit = Dpart;  % 200 X 496
% 
% thres=0.001;
% s=Dinit*trn;% added this line
% st=@(s,threshold)sign(s).*max(0,abs(s)-threshold);%added this line
% sparsed=st(s,thres);% added this line
% 
% W=HMat*sparsed' *inv(sparsed*sparsed' + eye(size(sparsed*sparsed')));
% dict.W = W;
% dict.D = Dinit;

% BATCH LEARNING
% M=496;
% theta = zeros(496,200);
% gamma = zeros(496,496);
% beta=0;
% bsize = num_images / 5;
% n0 = num_images/(bsize*10);
% for iter = 1 : num_iters
%     ind_rnd = randperm(num_images);
%     for batch = 1:(num_images/bsize)        
%         % load the dataset        
%         batch_idx = ind_rnd((1:bsize) + bsize*(batch-1));
%         yt = trn(:,batch_idx);
%         ht = HMat(:,batch_idx);
%         
%         thres=0.001;        
%         s=D*yt;% added this line
%         st=@(s,threshold)sign(s).*max(0,abs(s)-threshold);%added this line
%         sparsed=st(s,thres);
%         theta = (0.97)*theta + 0.3*(yt*sparsed');
% %         theta = theta./norm(theta);
%         gamma = (0.97)*gamma + 0.3*(yt*yt');
% %         gamma = gamma./norm(gamma);
%         beta = (0.97)*beta + 0.3*0.02*norm(yt,'fro');        
%         [U,S,V] = svd(gamma+beta*1);
%         Lj = U*sqrt(inv(S))*V;
%         Lj = Lj./norm(Lj);
%         [Q,S,R] = svd(Lj*theta);
%         D = (0.5*R')*(S+sqrt(S.^2+(2*beta*1)))'*(Q'*Lj);
%         dict.D = D;
%         
%         %W=HMat*sparsed' *inv(sparsed*sparsed' + eye(size(sparsed*sparsed')));
%         grad_W = (1-mu)*(W*sparsed - ht)*sparsed' + nu2*W;
%         rho_i = min(rho, rho*n0/batch);
%         W = W - rho_i*grad_W;
%         dict.W = W;
%     end
% end


% crf iterations
% for iter = 1 : num_iters
%     tic;
%     stat.fobj_total = 0;
%     % Take a random permutation of the samples
%     ind_rnd = randperm(num_images);
%     
%     for batch = 1:(num_images/bsize)
%         % load the dataset
%         % we only loads one sample or a small batch at each iteration
%         batch_idx = ind_rnd((1:bsize) + bsize*(batch-1));
%         yt = trn(:,batch_idx);
%         ht = HMat(:,batch_idx);
%         qt = QMat(:,batch_idx);
%         
%         D = model.D;
%         W = model.W;
%         A = model.A;
%         
%         % sparse coding
%         %S = L1QP_FeatureSign_Set(yt, D, gamma, lambda);
%         S = mexLasso(yt,D,param);
%        
%         % compute the gradient of crf parameters       
%         grad_W = (1-mu)*(W*S - ht)*S' + nu2*W; %
% %         grad_A = mu*(A*S - qt)*S' + nu1*A;
%         grad_S1 = W'*(W*S - ht); % gradient w.r.t S for 0.5*||H-WS||_2^2
%         grad_S2 = A'*(A*S - qt); % gradient w.r.t S for 0.5*||Q-AS||_2^2
%         
%         % compute the gradient of dictionary
%         % find the active set and compute beta
%         B1 = zeros(num_bases, bsize);
%         B2 = zeros(num_bases, bsize);
%         DtD = D'*D;
%         for j = 1:bsize
%             active_set = find(S(:,j) ~= 0);
%             %DtD = D(:,active_set)'*D(:,active_set) + gamma*eye(length(active_set));
%             DtD_hat = DtD(active_set,active_set) + gamma*eye(length(active_set));
%             
%             %DtD_inv = DtD\eye(length(active_set));
%             DtD_inv = DtD_hat\eye(length(active_set));
%             
%             beta1 = DtD_inv * grad_S1(active_set,j);
%             B1(active_set,j) = beta1;
%             
%             beta2 = DtD_inv * grad_S2(active_set,j);
%             B2(active_set,j) = beta2;
%         end
%         grad_D = (1-mu)*(-D*B1*S' + (yt - D*S)*B1') + mu*(-D*B2*S' + (yt - D*S)*B2'); % dD = -D*B*S' + (X - D*S)*B';
%         
%         %use yang's method
%         %gfullMat = zeros([size(D),size(D,2)]);
%         %[gMat, IDX] = sparseDerivative(D, full(S), yt);
%         %gfullMat(:,IDX,IDX) = gMat;
%         %gradSmat = repmat(reshape(grad_S1,[1 1 length(grad_S1)]),size(D));
%         %grad_D = sum(gfullMat.*gradSmat,3);
%         
%         % update the learning rate
%         rho_i = min(rho, rho*n0/batch);
%         
%         % update model parameters
%         D = D - rho_i*grad_D;
%         l2norm = (sum(D.^2)).^.5;
%         D = D ./ repmat(l2norm,size(D,1),1);
%         model.D = D;
%         
%         W = W - rho_i*grad_W;
%         model.W = W;
%         
% %         A = A - rho_i*grad_A;
% %         model.A = A;    
%     end
%   
%     % get statistics
%     S = mexLasso(trn,D,param);
%     fobj = get_objective(D, S, trn, W, HMat, A, QMat, lambda, mu);
%     stat.fobj_avg(iter) = fobj + 0.5*nu1*sum(sum(W.^2)) + 0.5*nu2*sum(sum(A.^2));
% %     fprintf('Iter = %d, Elapsed Time = %f\n', iter, toc);
% end

% figure(3);
% plot(stat.fobj_avg)