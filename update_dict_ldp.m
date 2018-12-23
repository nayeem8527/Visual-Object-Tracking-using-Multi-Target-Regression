function dict = update_dict_ldp(feat, dictpars, dict)

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



trn = feat.feaArr;
num_bases = dictpars.numBases;
num_iters = dictpars.maxIters;
num_images = size(trn, 2);
numClass = dictpars.numcls; % number of objects
HMat = zeros(numClass, num_bases);
HMat(1, feat.label == 1) = 1;
HMat(2, feat.label == -1) = 1;

%% compute Q
dictLabel = [];
numPerClass = num_bases / numClass;
for classid = 1:numClass
    labelvector = zeros(numClass,1);
    labelvector(classid) = 1;
    dictLabel = [dictLabel repmat(labelvector,1,numPerClass)];
end
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
nu2 = dictpars.nu2;
mu = dictpars.mu;
rho = dictpars.rho;
% bsize = dictpars.batchSize;
% bsize = num_images / 5;
% n0 = num_images/10
% n0 = num_images/(bsize*10);

Dpart = dict.D; % dictionary
W = dict.W; % classifier
% model.A = dict.A; % transform matrix

% param.lambda = dictpars.lambda;
% param.lambda2 = 0;
% param.mode = 2;

% DOING THE SAME WAY AS DONE IN INITIALIZATION
T0=11;  %sparsity level for each patch
l2=0.5; %weight for the negative log-determinant penalty in problem formulation
l4=l2;  %weight for the incoherence penalty in problem formulation
l3=l2;  %weight for the Frobenius norm penalty used within CG (in transform update step)
p=1;  %parameter for the incoherence penalty in problem formulation
mu=1e-9; %step size within CG in the transform update step
numt=5;  %number of iterations of alternating transform learning algorithm (i.e., iterations of sparse coding and transform update)
numg=10; %number of iterations of CG in the transform update step
Data = trn;
STY = T0*ones(1,size(Data,2));
[K,n]=size(Dpart); 
ix=find(STY>0); q=Data(:,ix); STY=STY(:,ix); N=size(q,2);
ez=K*(0:(N-1));STY=STY + ez;
ZX=Data*Data';
for i=1:numt    
    %%%%%%Sparse Coding Step%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    X1=Dpart*q;
    [s]=sort(abs(X1),'descend');
    X = X1.*(bsxfun(@ge,abs(X1),s(STY)));
    %%%%%%Transform Update Step%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    RSu=(q*X')';
    %CG iterations
    for j=1:numg
        ZZ=(Dpart*Dpart').^(p-1);
        cc = (-2*RSu) + (2*Dpart*ZX');
        cc2=(-2)*Dpart*(inv(Dpart'*Dpart));
        cc3= 2*Dpart;
        cc4=2*p*((ZZ*Dpart) - (((diag(ZZ))*(ones(1,n))).*Dpart));
        deD= l2*cc2 + cc + l3*cc3  + l4*cc4;
        if(j==1)
            g=deD;
            d=-g;
        else
            be = ((norm(deD,'fro'))^2)/((norm(g,'fro'))^2);
            g=deD;
            d= - g + be*d;
        end
        Dpart = Dpart + (mu)*d;   %constant step size
    end
    %post-normalization of the rows of the transform
    SL=diag(1./sqrt(sum((Dpart').^2)));
    Dpart=SL*Dpart;
end
% param1.D = Dinit;
Dinit = Dpart;  % 200 X 496
thres=0.001;
s=Dinit*trn;% added this line
st=@(s,threshold)sign(s).*max(0,abs(s)-threshold);%added this line
sparsed=st(s,thres);% added this line
W=HMat*sparsed' *inv(sparsed*sparsed' + eye(size(sparsed*sparsed')));
dict.W = W;
dict.D = Dinit;
