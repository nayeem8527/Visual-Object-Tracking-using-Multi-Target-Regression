function [bb_next, samples, pred,ind,gt] = do_tracking(img, dict,bb_prev, featpars, trackpars, sparsity, opt,p1,p2,gt)
%

%% global features
[bb_next, samples, pred,ind,gt] = use_global_feature(img, dict,bb_prev, trackpars, sparsity, opt,p1,p2,gt,featpars.bbox);