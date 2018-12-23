function feat = extract_global_feature(img, bbox, pos_num, neg_num,p1,p2,n1,n2,check)
%
% p = [(bbox(2)+bbox(4))/2, (bbox(1)+bbox(3))/2, bbox(4)-bbox(2), bbox(3)-bbox(1), 0];
% psize =  trackpars.nsize;
% psize = [32 32];
% p0 = [p(1), p(2), p(3)/psize(1), p(5), p(4)/p(3), 0]'; %param0 = [px, py, sc, th,ratio,phi];   
% p0 = affparam2mat(p0); 


% theta = [1,1,0,0,0,0];
for i=1:size(img,3)
    gfrm(:,:,i) = double(img(:,:,i));
end
% gfrm = double(img);

%%--- postive sample ---%%
% centers = repmat(affparam2geom(p0), [1, pos_num]);
% locs = centers + randn(6,pos_num) .* repmat(theta(:), [1, pos_num]);

% pos_examples = gen_samples('gaussian', bbox, 1, p1, p2,gfrm);

% r = overlap_ratio(pos_examples,bbox);
% pos_examples = pos_examples(r>0.7,:);
% pos_examples = pos_examples(randsample(end,min(pos_num/2,end)),:);

% wimgs = [];
% for i=1:size(pos_examples,1)
%     wimgs(:,:,i) = imcrop(gfrm,pos_examples(i,:));
% end
% wimgs = warpimg(gfrm, affparam2mat(locs), psize);                  

% feat.feaArr = [];
% feat.label = [ones(pos_num, 1); -1*ones(neg_num, 1)];

% p2 = gcp();
% F2 = parfeval(p2,@feat_extractor_hog,1,pos_num, wimgs);
feats_pos = feat_extractor_hog(gfrm,bbox);

%%--- negtive sample ---%%
% width= bbox(4) - bbox(2);
% height = bbox(3) - bbox(1);
% overlap_ratio = .05;
% cx = [width*(1-overlap_ratio) width*(1+overlap_ratio)];
% cy = [height*(1-overlap_ratio) height*(1+overlap_ratio)];
% centers(1,:) = centers(1,:) + random('uniform',cx(1),cx(2),1,neg_num) .* sign(rand(1,neg_num)-0.5);
% centers(2,:) = centers(2,:) + random('uniform',cy(1),cy(2),1,neg_num) .* sign(rand(1,neg_num)-0.5);
% wimgs_neg = warpimg(gfrm, affparam2mat(centers), psize);
% if check==1
%     neg_examples = [gen_samples('uniform', bbox, neg_num/2, n1, n2,gfrm);gen_samples('whole', bbox, neg_num/2, n1,n2,gfrm)];    
% %     r = overlap_ratio(neg_examples,bbox);
% %     neg_examples = neg_examples(r<0.2,:);
% %     neg_examples = neg_examples(randsample(end,min(neg_num,end)),:);
% else
%     neg_examples = [gen_samples('uniform', bbox, neg_num, n1, n2,gfrm)];
% %     r = overlap_ratio(neg_examples,bbox);
% %     neg_examples = neg_examples(r<0.2,:);
% %     neg_examples = neg_examples(randsample(end,min(neg_num/2,end)),:);
% end
% 
% feat.feaArr = [];
% feat.label = [ones(size(pos_examples,1), 1); -1*ones(size(neg_examples,1), 1)];
% r = overlap_ratio(neg_examples,bbox);
% neg_examples = neg_examples(r<0.5,:);
% neg_examples = neg_examples(randsample(end,min(300,end)),:);

% wimgs_neg = [];
% for i=1:size(pos_examples,1)
%     wimgs_neg(:,:,i) = imcrop(gfrm,neg_examples(i,:));
% end
% p3 = gcp();
% F3 = parfeval(p3,@feat_extractor_hog,1,neg_num, wimgs_neg);
% tic
% feats_neg = feat_extractor_hog(gfrm,neg_examples);
% val = fetchOutputs(F2);
% val1 = fetchOutputs(F3);
% feats_pos = F2.OutputArguments{1};
% feats_neg = F3.OutputArguments{1};            
% cancel(F2);
% cancel(F3);
feat.feaArr = [feats_pos];% feats_neg];
feat.boxes = [bbox];
% toc
A_norm = sqrt(sum(feat.feaArr .* feat.feaArr));
feat.feaArr = feat.feaArr ./ (ones(size(feat.feaArr,1),1) * A_norm + eps);
