function [bb_next, samples, pred,ind,gt] = use_global_feature(img, dict, bb_prev, trackpars, sparsity, opt,p1,p2,gt,ppp)
%
gfrm = double(img);

bbox = bb_prev;
% p = [(bbox(2)+bbox(4))/2, (bbox(1)+bbox(3))/2, bbox(4)-bbox(2), bbox(3)-bbox(1), 0];
% psize = trackpars.nsize;
% p0 = [p(1), p(2), p(3)/psize(1), p(5), p(4)/p(3), 0]'; %param0 = [px, py, sc, th,ratio,phi];   
% p0 = affparam2mat(p0); 

%% generate candidates
% pos_examples = gen_samples('gaussian', bbox, 500, 0.3, 2,gfrm);
pos_examples = gen_samples('gaussian', bbox, 1000, p1, p2,gfrm);
% r = overlap_ratio(pos_examples,bbox);
% pos_examples = pos_examples(r>0.7,:);
% pos_examples = pos_examples(pos_examples(:,1)>=1,:);
% pos_examples = pos_examples((pos_examples(:,1)+pos_examples(:,3))<=size(img,2),:);
% pos_examples = pos_examples(pos_examples(:,2)>=1,:);
% pos_examples = pos_examples((pos_examples(:,2)+pos_examples(:,4))<=size(img,1),:);
% pos_examples = pos_examples(randsample(end,min(1000,end)),:);

% neg_examples = gen_samples('uniform', bbox, 500, 1, 5,gfrm);
% %     gen_samples('whole', bbox, 200, 1,10,gfrm)];
% r = overlap_ratio(neg_examples,bbox);
% neg_examples = neg_examples(r<0.2,:);
% neg_examples = neg_examples(neg_examples(:,1)>=1,:);
% neg_examples = neg_examples((neg_examples(:,1)+neg_examples(:,3))<=size(img,2),:);
% neg_examples = neg_examples(neg_examples(:,2)>=1,:);
% neg_examples = neg_examples((neg_examples(:,2)+neg_examples(:,4))<=size(img,1),:);
% 
% neg_examples = neg_examples(randsample(end,min(500,end)),:);

examples = [pos_examples];%; neg_examples];
% locs = sampling(p0, opt.numsample, opt.affsig); 

% candidates = warpimg(gfrm, affparam2mat(locs), psize);
% candidates = [];
% for i=1:size(examples,1)
%     candidates(:,:,i) = imcrop(gfrm,examples(i,:));
% end

disp('use global feature');
% tic
% data = [];
% parfor i = 1:length(locs)
%     sample = candidates(:,:,i);
% 
%     %% hog features
% %     hog = double(vl_hog(im2single(sample), 8));
% %     feat.feaArr = [feat.feaArr double(hog)];
%     data = [data hog(:)];
% end
% 
% options.scale         = [1,2,3];
% options.deltax        = 12;
% options.deltay        = 12;
% options.n             = 2;
% options.ldporient     = [0 , 1 , 2 , 3];
% options.color         = 0;
% options.patchsize     = 20;
% options.rmextremebins = 1;
% options.norm          = 4;
% options.clamp         = 0.2;
% 
% data_ldp = [];
% parfor i = 1:length(locs)
%     sample = candidates(:,:,i);
%     [hog , ~] =  denseMBLDP(uint8(sample) , options );
%     hog = hog';
%     [coeff,~] = pca(hog);
%     coeff = coeff(:,1:300);
%     hog = hog*coeff;
%     hog = sum(hog,2);
%     hog = (hog-min(hog))/(max(hog)-min(hog));
%     data_ldp = [data_ldp hog(:)];
% end
% toc

% breakings = 1:ceil(length(locs)/4):length(locs);
data = feat_extractor_hog(gfrm,examples,ppp);
% p6 = gcp();
% F6 = parfeval(p6,@feat_extractor_hog,1,length(locs), candidates);
% p71 = gcp();
% F71 = parfeval(p71,@feat_extractor_ldp,1,breakings(1),breakings(2), candidates);
% p72 = gcp();
% F72 = parfeval(p72,@feat_extractor_ldp,1,breakings(2)+1,breakings(3), candidates);
% p73 = gcp();
% F73 = parfeval(p73,@feat_extractor_ldp,1,breakings(3)+1,breakings(4), candidates);
% p74 = gcp();
% F74 = parfeval(p74,@feat_extractor_ldp,1,breakings(4)+1,length(locs), candidates);
% value = fetchOutputs(F6);
% value1 = fetchOutputs(F71);
% value2 = fetchOutputs(F72);
% value3 = fetchOutputs(F73);
% value4 = fetchOutputs(F74);
% data = F6.OutputArguments{1};
% data_ldp_1 = F71.OutputArguments{1};            
% data_ldp_2 = F72.OutputArguments{1};            
% data_ldp_3 = F73.OutputArguments{1};            
% data_ldp_4 = F74.OutputArguments{1};            
% cancel(F6);
% cancel(F71);
% cancel(F72);
% cancel(F73);
% cancel(F74);
% toc

A_norm = sqrt(sum(data .* data));
data = data ./ (ones(size(data,1),1) * A_norm + eps);

% A_norm_ldp_1 = sqrt(sum(data_ldp_1 .* data_ldp_1));
% data_ldp_1 = data_ldp_1 ./ (ones(size(data_ldp_1,1),1) * A_norm_ldp_1 + eps);
% A_norm_ldp_2 = sqrt(sum(data_ldp_2 .* data_ldp_2));
% data_ldp_2 = data_ldp_2 ./ (ones(size(data_ldp_2,1),1) * A_norm_ldp_2 + eps);
% A_norm_ldp_3 = sqrt(sum(data_ldp_3 .* data_ldp_3));
% data_ldp_3 = data_ldp_3 ./ (ones(size(data_ldp_3,1),1) * A_norm_ldp_3 + eps);
% A_norm_ldp_4 = sqrt(sum(data_ldp_4 .* data_ldp_4));
% data_ldp_4 = data_ldp_4 ./ (ones(size(data_ldp_4,1),1) * A_norm_ldp_4 + eps);
% data_ldp = [data_ldp_1 data_ldp_2 data_ldp_3 data_ldp_4];


disp('Classify');
% tic
[pred,GT] = classify(double(data),dict, gt);
% toc
% [~, ind] = min(pred.val);
[val, ind] = max(pred);
disp(val);
gt = GT(:,ind);
% ind = pred.best_ind;
% best = locs(:, ind);
best = examples(ind,:);
% bbox = [best(1), best(2), best(3)*psize(2), best(5)*best(3)*psize(1)];

% bb_next = [bbox(2)-bbox(4)/2, bbox(1)-bbox(3)/2, bbox(2)+bbox(4)/2, bbox(1)+bbox(3)/2];
bb_next=best;

isgood = update_check(pred, ind, dict);

samples = [];
if isgood
    disp('again extract global features');
%     tic
%     samples = extract_global_feature(img, bb_next, trackpars.updatenum, trackpars.updatenum);
    samples = extract_global_feature(img, bb_next, 100,100,0.1,3,2,3,2,ppp);
%     samples_ldp = extract_global_feature_ldp(img, bb_next, trackpars.updatenum, trackpars.updatenum);
%     toc
end

