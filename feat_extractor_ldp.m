function feats = feat_extractor_ldp(start,num,wimgs)
   
options.scale         = [1,2,3];
options.deltax        = 12;
options.deltay        = 12;
options.n             = 2;
options.ldporient     = [0 , 1 , 2 , 3];
options.color         = 0;
options.patchsize     = 20;
options.rmextremebins = 1;
options.norm          = 4;
options.clamp         = 0.2;

feats = [];
for i = start:num
    sample = wimgs(:,:,i);
    [hog , ~] =  denseMBLDP(uint8(sample) , options );
    hog = hog';
    [coeff,~] = pca(hog);
    coeff = coeff(:,1:300);
    hog = hog*coeff;
    hog = sum(hog,2);
    hog = (hog-min(hog))/(max(hog)-min(hog));
    feats = [feats hog(:)];
end
end