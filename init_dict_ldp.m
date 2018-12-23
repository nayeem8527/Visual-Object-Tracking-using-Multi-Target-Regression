function [Dinit, Winit, Ainit, Q] = init_dict_ldp(feat, dictpars)
%

T0=11;  %sparsity level for each patch
l2=0.5; %weight for the negative log-determinant penalty in problem formulation
l4=l2;  %weight for the incoherence penalty in problem formulation
l3=l2;  %weight for the Frobenius norm penalty used within CG (in transform update step)
p=1;  %parameter for the incoherence penalty in problem formulation
mu=1e-9; %step size within CG in the transform update step
numt=5;  %number of iterations of alternating transform learning algorithm (i.e., iterations of sparse coding and transform update)
numg=10; %number of iterations of CG in the transform update step

trn = feat.feaArr;
dictsize = dictpars.numBases;
% iterations = dictpars.iterationini;
numClass = dictpars.numcls; % number of objects
Dinit = []; % for C-Ksvd and D-Ksvd
dictLabel = [];
numPerClass = dictsize / numClass;
H_train = zeros(numClass, dictsize);
H_train(1, feat.label == 1) = 1;
H_train(2, feat.label == -1) = 1;

for classid = 1:numClass
    col_ids = find(H_train(classid,:) == 1);
    data_ids = find(colnorms_squared_new(trn(:, col_ids)) > 1e-6);   % ensure no zero data elements are chosen
    perm = randperm(length(data_ids));
    
    %% Initilization for LC-KSVD (perform KSVD in each class)
    Dpart = trn(:, col_ids(data_ids(perm(1:min(numPerClass,length(perm))))))'; %initialized dictionary 100 X 496
    Data = trn(:,col_ids(data_ids));
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
%     param1.mode = 2;
%     param1.K = dictsize;
%     param1.lambda = dictpars.lambda;
%     param1.lambda2 = 0;
%     param1.iter = iterations;
%     param1.D = Dpart;
%     Dpart = mexTrainDL(trn(:,col_ids(data_ids)), param1);
    Dinit = [Dinit;Dpart];  % 200 X 496
    labelvector = zeros(numClass,1);
    labelvector(classid) = 1;
    dictLabel = [dictLabel repmat(labelvector,1,numPerClass)];
end

Dpart = Dinit;
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

% param2.lambda = dictpars.lambda;
% param2.lambda2 = 0;
% param2.mode = 2;
% Xinit = mexLasso(trn, Dinit, param2);
thres=0.001;
s=Dinit*trn;% added this line
st=@(s,threshold)sign(s).*max(0,abs(s)-threshold);%added this line
Xinit=st(s,thres);% added this line

% ISTA
% alpha = max(max(eig(Dinit'*Dinit)));
% T = 0.5/(2*alpha);
% x = 0*Dinit'*trn; % 200 x 400
% for k = 1:10    
%     Hx = Dinit*x;
%     x = wthresh(x + (Dinit'*(trn - Hx))/alpha,'s', T);
% end
% Xinit = x;
% learning linear classifier parameters
% Winit = inv(Xinit*Xinit' + eye(size(Xinit*Xinit')))*Xinit*H_train';
% Winit = Winit';
Winit=H_train*Xinit' *inv(Xinit*Xinit' + eye(size(Xinit*Xinit')));
Q = zeros(dictsize, size(trn, 2)); % energy matrix
for frameid = 1:size(trn, 2)
    label_training = H_train(:, frameid);
    [~, maxid1] = max(label_training);
    for itemid = 1:size(Dinit',2)
        label_item = dictLabel(:, itemid);
        [~, maxid2] = max(label_item);
        if(maxid1 == maxid2)
            Q(itemid, frameid) = 1;
        end
    end
end
% 
% Ainit = inv(Xinit*Xinit' + eye(size(Xinit*Xinit')))*Xinit*Q';
% Ainit = Ainit';
Ainit = eye(size(Xinit,1));