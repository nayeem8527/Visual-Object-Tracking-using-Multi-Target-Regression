function [dict, dictpars, bbox] = init_tracker(img, featpars, dictpars, trackpars)
%
bbox = featpars.bbox;
% bbox = [p(2)-p(4)/2 p(1)-p(3)/2 p(2)+p(4)/2 p(1)+p(3)/2];
% width= bbox(4) - bbox(2);
% height = bbox(3) - bbox(1);
% width = bbox(3);
% height = bbox(4);
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

% region = imresize(region, featpars.dims);
% gt = imresize(gt, featpars.dims);

%% extract features
pos_num = trackpars.posnum;
neg_num = trackpars.negnum;
disp('extracting the global features');
% feat = extract_global_feature(img, bbox, pos_num, neg_num);
% feat = extract_global_feature(img, bbox, 1000, 5000,0.1,5,1,10,1);
feat = extract_global_feature(img, bbox, 500, 500,0.1,5,1,10,1);

%% dictionary initialization
% numcls = length(unique(feat.label));
% dictpars.numBases = dictpars.numpercls * numcls; % dictionary size
dictpars.numcls = 1;
disp('initializing the dictionary');
[D, W] = init_dict(feat, dictpars);
dict.D = D;
dict.W = W;
% dict.A = A;
% dict.Q = Q;
bbox = bbox(10,:);
% T = []; TW = [];
% [T, TW] = get_template(img, bbox, T, TW, 0, trackpars);
% dict.T = T;
% dict.TW = TW;