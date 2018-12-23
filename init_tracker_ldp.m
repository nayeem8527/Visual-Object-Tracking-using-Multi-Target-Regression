function [dict_ldp,dictpars,bbox] = init_tracker_ldp(img, featpars, dictpars, trackpars)
%
bbox = featpars.bbox;
% bbox = [p(2)-p(4)/2 p(1)-p(3)/2 p(2)+p(4)/2 p(1)+p(3)/2];
% width= bbox(4) - bbox(2);
% height = bbox(3) - bbox(1);
% topx = max(1, bbox(2) - width/2);
% topy = max(1, bbox(1) - height/2);
% bottomx = min(size(img,2), bbox(4) + width/2);
% bottomy = min(size(img,1), bbox(3) + height/2);
% 
% gt = zeros(size(img,1), size(img,2));
% gt(bbox(1):bbox(3), bbox(2):bbox(4)) = 1;
% 
% region = img(topy:bottomy, topx:bottomx, :);
% gt = gt(topy:bottomy, topx:bottomx);
% 
% region = imresize(region, featpars.dims);
% gt = imresize(gt, featpars.dims);

%% extract features
pos_num = trackpars.posnum;
neg_num = trackpars.negnum;
disp('extracting the global features');
feat = extract_global_feature_ldp(img, bbox, pos_num, neg_num);

%% dictionary initialization
numcls = length(unique(feat.label));
dictpars.numBases = dictpars.numpercls * numcls; % dictionary size
dictpars.numcls = numcls;
disp('initializing the dictionary');
[D, W, A, Q] = init_dict_ldp(feat, dictpars);
dict_ldp.D = D;
dict_ldp.W = W;
dict_ldp.A = A;
dict_ldp.Q = Q;

T = []; TW = [];
[T, TW] = get_template_ldp(img, bbox, T, TW, 0, trackpars);
dict_ldp.T = T;
dict_ldp.TW = TW;