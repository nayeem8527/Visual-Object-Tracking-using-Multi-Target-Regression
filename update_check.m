function isgood = update_check(pred, ind, dict)
%
isgood = 1;
% numcls = 2;
% th_cls = 0.65;
% th_sc = 0.3;
% numperclass = size(dict.D,2) / numcls;
% sc = pred.sc(:,ind);

% for i = 1:numcls
%     prob(i) = sum(abs(sc(numperclass*(i-1)+1:i*numperclass))) / sum(abs(sc));
% end
% 
% ent = -sum(prob .* log(prob));


% if err > th_sc || ent > th_cls 
%     isgood = 0;
% end
% 

% avgW = sum(dict.TW) / length(dict.TW);
% th_lowerbound = .95;

% pred.cls_err(ind)
% pred.rec_err(ind) 
% disp('Reconstruction error:')
% disp(pred.rec_err(ind));
% disp(pred.rec_err_ldp(ind));
% disp('Classification error:')
% disp(pred.cls_err(ind));
% disp(pred.rec_err(ind)+pred.cls_err(ind));
% disp(pred.cls_err_ldp(ind));
% if max(pred)<0.55  %pred.cls_err(ind) > 0.6 || pred.rec_err(ind) > 0.45
%     isgood = 0;
% end
% if exp(-pred.val(ind)) < th_lowerbound * avgW
%     isgood = 0;
% end