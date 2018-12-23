clear all; close all;
clc
warning('off');

% cd ./sampling
% mex -O interp2.cpp
% cd ..
    
addpath('sampling');
% addpath('model','matconvnet/matlab');
% 
% vl_setupnn();
% vl_compilenn('enableGpu',false,'EnableImreadJpeg',false);
disp('libraries compiled');

%% change the directory of test sequences
seq_path = 'data/';

s = dir(seq_path);
seq = s(3:end);
%% loop through all sequences
for n = 1:length(seq)
    trackpars.title = seq(n).name;
    
    %% parameter setting
    trackparam;

    %% collect training samples
    dataPath = strcat(dataPath,'/');
    img = imread([dataPath files(10).name]);
    if size(img,3) == 3
        img = rgb2gray(img);
    end
    
    
     for i=1:10
         img = imread([dataPath files(i).name]);
         img = rgb2gray(img);
         I(:,:,i) = img;
     end
    %%?initialize tracker 
    disp('initializing the tracker');
%     tic    
    [dict, dictpars, bb_prev] = init_tracker(I, featpars, dictpars, trackpars);                            
%     [dict_ldp, ~, ~] = init_tracker_ldp(img, featpars, dictpars, trackpars);                       
%     toc
%     tic
%     p = gcp();
%     F = parfeval(p,@init_tracker,3,img,featpars,dictpars,trackpars);
%     p1 = gcp();
%     F1 = parfeval(p1,@init_tracker_ldp,3,img,featpars,dictpars,trackpars);
%     value = fetchOutputs(F);
%     value1 = fetchOutputs(F1);
%     dict = F.OutputArguments{1};
%     dictpars = F.OutputArguments{2};
%     bb_prev = F.OutputArguments{3};
%     dict_ldp = F1.OutputArguments{1};
%     cancel(F);
%     cancel(F1);
%     toc
    
    drawopt = [];
    samples.feaArr = [];
%     samples.label = [];
    samples.boxes = [];
%     samples_ldp.feaArr = [];
%     samples_ldp.label = [];
    result = bb_prev;
    cnt = 0;

    %% draw results

    drawopt = draw_result(drawopt, 1, img, bb_prev); 
    imwrite(frame2im(getframe(gcf)),sprintf('results/%s/%s_global_0001.png', trackpars.title,  trackpars.title));


    tic;
    duration = 0;
    %% do tracking 
    gt = [1;0];
    p1=0.6;
%     p1=0.3;
    for f = 11 : frameNum
%         try
            imgName = sprintf('%s08d.jpg', dataPath, f); 
            if exist(imgName, 'file')
                frame = imread(imgName);
            else
                imgName = sprintf('%s%08d.jpg', dataPath, f); 
                frame = imread(imgName);
            end

            if size(frame,3) == 3
                frame = rgb2gray(frame);
            end
            disp('performing the tracking'); 
%             if f==2                
            [bb_next, feats,gt] = do_tracking(frame, dict, bb_prev, featpars, trackpars, dictpars.sparsity, opt,p1,1,gt);
%             else
%                 if (pred.rec_err(ind)+pred.cls_err(ind))>1                    
%                      p1 = min(1.5, 1.1*p1);                     
%                      [bb_next, feats, pred,ind] = do_tracking(frame, dict, bb_prev, featpars, trackpars, dictpars.sparsity, opt,p1,1);        
%                 else
%                     [bb_next, feats, pred,ind] = do_tracking(frame, dict, bb_prev, featpars, trackpars, dictpars.sparsity, opt,p1,1);        
%                 end                
%             end

            if ~isempty(feats)
                cnt = cnt + 1;
                samples.feaArr = [samples.feaArr feats];
%                 samples.label = [samples.label feats.label'];
                samples.boxes = [samples.boxes;gt];
    %             samples_ldp.feaArr = [samples_ldp.feaArr feats_ldp.feaArr];
    %             samples_ldp.label = [samples_ldp.label feats_ldp.label'];
%                 disp('get template');
%                 tic
%                 [dict.T, dict.TW] = get_template(frame, bb_next, dict.T, dict.TW, max(pred), trackpars);
    %             [dict_ldp.T, dict_ldp.TW] = get_template_ldp(frame, bb_next, dict_ldp.T, dict_ldp.TW, pred.t2(ind), trackpars);
    %             toc
%                 p10 = gcp();
%                 F10 = parfeval(p10,@get_template,2,frame, bb_next, dict.T, dict.TW, pred.val(ind), trackpars);
    %             p11 = gcp();
    %             F11 = parfeval(p11,@get_template_ldp,2,frame, bb_next, dict_ldp.T, dict_ldp.TW, pred.t2(ind), trackpars);
%                 value = fetchOutputs(F10);
    %             value1 = fetchOutputs(F11);
%                 dict.T = F10.OutputArguments{1};
%                 dict.TW = F10.OutputArguments{2};
    %             dict_ldp.T = F11.OutputArguments{1};
    %             dict_ldp.TW = F11.OutputArguments{2};
%                 cancel(F10);
    %             cancel(F11);
%                 toc
            end

            if ~isempty(samples.feaArr) && size(samples.feaArr,2)==5 && trackpars.isupdate ~= 0            
                
                dict = update_dict(samples, dictpars, dict);
                
                
    %             dict_ldp = update_dict_ldp(samples_ldp, dictpars, dict_ldp);
%                 disp('update dict');
%                 tic
%                 p12 = gcp();
%                 F12 = parfeval(p12,@update_dict,1,samples, dictpars, dict);
    %             p13 = gcp();
    %             F13 = parfeval(p13,@update_dict_ldp,1,samples_ldp, dictpars, dict_ldp);
%                 value = fetchOutputs(F12);
    %             value1 = fetchOutputs(F13);
%                 dict = F12.OutputArguments{1};
    %             dict_ldp = F13.OutputArguments{1};            
%                 cancel(F12);
    %             cancel(F13);
%                 toc
                samples.feaArr = [];
%                 samples.label = [];   
                samples.boxes = [];
    %             samples_ldp.feaArr = [];
    %             samples_ldp.label = [];            
            end

            %% draw result        
            drawopt = draw_result(drawopt, f, frame, bb_next); 
            imwrite(frame2im(getframe(gcf)),sprintf('results/%s/%s_global_%04d.png', trackpars.title, trackpars.title, f));   

            %% append new result
            result = [result; bb_next]; 
            bb_prev = bb_next;        
%         catch
%             disp('Got Error');
%             break;
%         end
    end            

    duration = duration + toc;      
    fprintf('%d frames took %.3f seconds : %.3fps\n', f, duration, f/duration);
    fps = f/duration;

    %% save results to file
    fileName = sprintf('results/%s_results.txt', trackpars.title);
    dlmwrite(fileName, result, 'delimiter', ' ');
    close all
 end