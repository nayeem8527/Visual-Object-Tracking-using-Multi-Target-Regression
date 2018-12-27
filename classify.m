function [pred,GT] = classify(data,dict,gt)

% HOG
% p8 = gcp();
% F8 = parfeval(p8,@classification,2,data,dict);

% % LDP
% p9 = gcp();
% F9 = parfeval(p9,@classification,2,data_ldp,dict_ldp);

% value = fetchOutputs(F8);
% value1 = fetchOutputs(F9);
% rec_err = F8.OutputArguments{1};
% cls_err = F8.OutputArguments{2};            
% rec_err_ldp = F9.OutputArguments{1};
% cls_err_ldp = F9.OutputArguments{2};            
% cancel(F8);
% cancel(F9);

% [rec_err,cls_err] = classification(data,dict);
[pred,GT] = classification(data,dict,gt);
% [rec_err_ldp,cls_err_ldp] = classification(data_ldp,dict_ldp);


% w = 0.5*(1-(rec_err./(rec_err+rec_err_ldp)));
% w_ldp = 0.5*(1-(rec_err_ldp./(rec_err+rec_err_ldp)));

% t = dict.W * Gamma;
% cls_err = exp(-(t(1,:)-t(2,:)));
% pred = [];
% pred.rec_err = rec_err;
% pred.cls_err = cls_err;

% pred.rec_err_ldp = rec_err_ldp;
% pred.cls_err_ldp = cls_err_ldp;

% % combine
% t1 = rec_err + w.*cls_err;
% pred.t1 = t1;
% t2 = rec_err_ldp + w_ldp.*cls_err_ldp;
% pred.t2 = t2;
% alpha = 0.8;
% pred.val = alpha * rec_err/sum(rec_err) + (1-alpha) * cls_err/sum(cls_err);
% pred.val = t1+t2;
