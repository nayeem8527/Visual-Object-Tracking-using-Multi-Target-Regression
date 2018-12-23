function feat = extract_global_feature_ldp(img, bbox, pos_num, neg_num)
%
p = [(bbox(2)+bbox(4))/2, (bbox(1)+bbox(3))/2, bbox(4)-bbox(2), bbox(3)-bbox(1), 0];
% psize =  trackpars.nsize;
psize = [32 32];
p0 = [p(1), p(2), p(3)/psize(1), p(5), p(4)/p(3), 0]'; %param0 = [px, py, sc, th,ratio,phi];   
p0 = affparam2mat(p0); 


theta = [1,1,0,0,0,0];
gfrm = double(img);

%%--- postive sample ---%%
centers = repmat(affparam2geom(p0), [1, pos_num]);
locs = centers + randn(6,pos_num) .* repmat(theta(:), [1, pos_num]);
wimgs = warpimg(gfrm, affparam2mat(locs), psize);                  
feat.feaArr = [];
feat.label = [ones(pos_num, 1); -1*ones(neg_num, 1)];

p4 = gcp();
F4 = parfeval(p4,@feat_extractor_ldp,1,1,pos_num, wimgs);
% feats_pos = feat_extractor_ldp(pos_num,wimgs);
% toc

%%--- negtive sample ---%%
width= bbox(4) - bbox(2);
height = bbox(3) - bbox(1);
overlap_ratio = .05;
cx = [width*(1-overlap_ratio) width*(1+overlap_ratio)];
cy = [height*(1-overlap_ratio) height*(1+overlap_ratio)];
centers(1,:) = centers(1,:) + random('uniform',cx(1),cx(2),1,neg_num) .* sign(rand(1,neg_num)-0.5);
centers(2,:) = centers(2,:) + random('uniform',cy(1),cy(2),1,neg_num) .* sign(rand(1,neg_num)-0.5);
wimgs = warpimg(gfrm, affparam2mat(centers), psize);

% feats_neg = feat_extractor_ldp(neg_num,wimgs);
p5 = gcp();
F5 = parfeval(p5,@feat_extractor_ldp,1,1,neg_num, wimgs);

value = fetchOutputs(F4);
value1 = fetchOutputs(F5);
feats_pos = F4.OutputArguments{1};
feats_neg = F5.OutputArguments{1};            
cancel(F4);
cancel(F5);
feat.feaArr = [feats_pos feats_neg];

A_norm = sqrt(sum(feat.feaArr .* feat.feaArr));
feat.feaArr = feat.feaArr ./ (ones(size(feat.feaArr,1),1) * A_norm + eps);
