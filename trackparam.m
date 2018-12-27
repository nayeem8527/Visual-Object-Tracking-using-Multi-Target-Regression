
%% parameters for extracting features
featpars.gridSpacing = 4;
featpars.patchSizes = [16];
featpars.patchSize = 16;
featpars.numOBins = 16;
featpars.numSBins = 4;
featpars.maxImSize = 640;
featpars.dims = [128 128];

%% parameters of dictionary learning
dictpars.gamma = 1e-6;
dictpars.lambda = 0.5;
dictpars.mu = 0.6; % ||Q-AX||^2
dictpars.nu1 = 1e-6; % regularization of A
dictpars.nu2 = 1e-6; % regularization of W
dictpars.rho = .2; % initial learning rate
dictpars.maxIters = 30; % iteration number for incremental learning
dictpars.batchSize = 100;
dictpars.iterationini = 5; % iteration number for initialization
dictpars.numpercls = 100;
dictpars.sparsity = 10;

%% tracking parameters 
trackpars.isupdate = 1;
trackpars.update = 5;
trackpars.posnum = 200;
trackpars.negnum = 200;
trackpars.updatenum = 100;
trackpars.nsize = [32 32];
trackpars.lengthT = 20;


%% individual parameters
switch (trackpars.title)
    case 'fsfsd';
        disp('Go to otherwise');
  
%     case 'Basketball';
%         opt = struct('numsample',800, 'affsig',[4,4,.08,.000,.0005,.000]);
%         featpars.bbox =[198+17 214+40 34 81];
%     
%     case 'BlurBody';
%         opt = struct('numsample',800, 'affsig',[4,4,.08,.000,.0005,.000]);
%         featpars.bbox =[400+87/2 48+319/2 87 319];
%         
%     case 'Car1';
%         opt = struct('numsample',800,'affsig',[5,5,.05,.00,.001,.001]);
%         featpars.bbox =[23+33 88+55/2 66 55];
%         
%     case 'Soccer';
%         opt = struct('numsample',800, 'affsig',[4,4,.08,.000,.0005,.000]);
%         featpars.bbox =[302+67/2 135+81/2 67 81];%Soccer
%     
%     case 'Singer2';
%         opt = struct('numsample',800, 'affsig',[4,4,.08,.000,.0005,.000]);
%         featpars.bbox =[298+67/2 149+122/2 67 122];
%                
%     case 'CarDark';
%         opt = struct('numsample',800,'affsig',[5,5,.05,.00,.001,.001]);
%         featpars.bbox =[73+29/2 126+23/2 29 23];
%         
%     case 'Coke';
%         opt = struct('numsample',800, 'affsig',[3,3, 0.03, 0.00, 0.0, .00]);
%         featpars.bbox =[298+24 160+40 48 80]  ;
%     
%     case 'Couple';
%         opt = struct('numsample',800, 'affsig',[3,3, 0.03, 0.00, 0.0, .00]);
%         featpars.bbox =[51+25/2 47+62/2  25 62];
%         
%     case 'Crowds';
%         opt = struct('numsample',800, 'affsig',[3,3, 0.03, 0.00, 0.0, .00]);
%         featpars.bbox =[561+11 311+51/2 22 51];
%     
%     case 'David2';
%         opt = struct('numsample',800, 'affsig',[7,7,.015,.000,.001,.001]);
%         featpars.bbox =[141+27/2 73+17 27 34];%David2
%         
%     case 'Deer';
%         opt = struct('numsample',800, 'affsig',[22 ,22, 0.03, 0.00, 0.0, .00]);
%         featpars.bbox = [350 40 100 70];
%         
%     case 'Dog1';
%         opt = struct('numsample',800, 'affsig',[22 ,22, 0.03, 0.00, 0.0, .00]);
%         featpars.bbox =[139+51/2 112+36/2 51 36];%Dog1
%         
%     case 'Fish';
%         opt = struct('numsample',800, 'affsig',[22 ,22, 0.03, 0.00, 0.0, .00]);
%         featpars.bbox =[134+60/2 55+88/2 60 88];
%         
%     case 'Ironman';
%         opt = struct('numsample',800, 'affsig',[22 ,22, 0.03, 0.00, 0.0, .00]);
%         featpars.bbox =[206+49/2 85+57/2 49 57];%Ironman
%         
%     case 'Jumping';
%         opt = struct('numsample',800, 'affsig',[22 ,22, 0.03, 0.00, 0.0, .00]);
%         featpars.bbox =[147+17 110+33/2 34 33];%Jumping
%         
%     case 'Liquor';
%         opt = struct('numsample',800, 'affsig',[22 ,22, 0.03, 0.00, 0.0, .00]);
%         featpars.bbox =[256+73/2 152+210/2 73 210];
%         
%     case 'Matrix';
%         opt = struct('numsample',800, 'affsig',[22 ,22, 0.03, 0.00, 0.0, .00]);
%         featpars.bbox =[331+38/2 39+42/2 38 42];%Matrix
%         
%     case 'MountainBike';
%         opt = struct('numsample',800, 'affsig',[22 ,22, 0.03, 0.00, 0.0, .00]);
%         featpars.bbox =[319+67/2 185+56/2 67 56];%Mountain Bike
%         
%     case 'Skating1';
%         opt = struct('numsample',800, 'affsig',[4,4,.08,.000,.0005,.000]);
%         featpars.bbox =[162+34/2 188+84/2 34 84];%Skating1
%         
%     case 'Subway';
%         opt = struct('numsample',800, 'affsig',[4,4,.08,.000,.0005,.000]);
%         featpars.bbox =[16+19/2 88+51/2 19 51];
%         
%     case 'Tiger1';
%         opt = struct('numsample',800, 'affsig',[4,4,.08,.000,.0005,.000]);
%         featpars.bbox =[232+76/2 88+84/2 76 84];%Tiger2
%         
%     case 'singer1';
%         opt = struct('numsample',800, 'affsig',[4,4,.08,.000,.0005,.000]);
%         featpars.bbox = [100 200 70 270];
%    
%     case 'davidin300';
%         opt = struct('numsample',800, 'affsig',[7,7,.015,.000,.001,.001]);    
%         featpars.bbox = [160 112 60 92];
%               
%     case 'car4';  
%         opt = struct('numsample',800,'affsig',[5,5,.05,.00,.001,.001]);                          
%         featpars.bbox = [245 180 200 150];
%         
%     case 'stone';
%         opt = struct('numsample',800, 'affsig',[6,6,.03,eps,eps,eps]);
%         featpars.bbox = [115 150 43 20];    
% 
%     case 'Bolt';
%         opt = struct('numsample',800, 'affsig',[6,5,.001,.001,.001,.000]);
%         featpars.bbox = [350 195 25 60];
%         
%     case 'girl'; 
%         opt = struct('numsample',800,'affsig',[6,4, 0.013, 0.00, 0.001, .00]);
%         featpars.bbox = [319 204 38 150];
%                 
%     case 'animal';
%         opt = struct('numsample',800, 'affsig',[22 ,22, 0.03, 0.00, 0.0, .00]);
%         featpars.bbox = [350 40 100 70];
% 
%     case 'Football';
%         opt = struct('numsample',800, 'affsig',[4,4, 0.02, 0.00, 0.0, .00]);
%         featpars.bbox = [330 125 50 50];
%    
%     case 'twinnings';
%         opt = struct('numsample',800, 'affsig',[3,3, 0.03, 0.00, 0.0, .00]);
%         featpars.bbox = [163,192,73,53];
%     
%     case 'Biker';
%         opt = struct('numsample',800, 'affsig',[4,4,.08,.000,.0005,.000]);
%         featpars.bbox =[262+16/2 94+13 16 26];
%     
%     case 'BlurCar2';
%         opt = struct('numsample',800, 'affsig',[4,4,.08,.000,.0005,.000]);
%         featpars.bbox =[227+122/2 207+99/2 122 99];
%         
%     % vot 2017
%     case 'ants1';
%         opt = struct('numsample',800,'affsig',[4,4,.08,.000,.0005,.000]);
%         dataPath = [seq_path trackpars.title '/'];
%         gt = importdata([dataPath 'groundtruth.txt']);
%         x = gt(:,1:2:end);
%         y = gt(:,2:2:end);
%         gt = [min(x,[],2), min(y,[],2), max(x,[],2) - min(x,[],2), max(y,[],2) - min(y,[],2)];
%         gt = gt(1,:);
% %         gt(1,1) = gt(1,1)+(gt(1,3)/2);
% %         gt(1,2) = gt(1,2)+(gt(1,4)/2);
%         featpars.bbox = gt;            
           
    otherwise;
        opt = struct('numsample',800,'affsig',[4,4,.08,.000,.0005,.000]);
        dataPath = [seq_path trackpars.title '/'];
        gt = importdata([dataPath 'groundtruth_rect.txt']);
%         x = gt(:,1:2:end);
%         y = gt(:,2:2:end);
%         gt = [min(x,[],2), min(y,[],2), max(x,[],2) - min(x,[],2), max(y,[],2) - min(y,[],2)];
        gt = gt(1,:);
%         gt(1,1) = gt(1,1)+(gt(1,3)/2);
%         gt(1,2) = gt(1,2)+(gt(1,4)/2);
        featpars.bbox = gt;            
%         error(['unknown title ' trackpars.title]);
end

%% path of the sequences, you can change it to the directory that contains the sequences 
dataPath = [seq_path trackpars.title '/'];

%% get frame number and resolution
files = dir([dataPath '/img/*.jpg']);
frameNum = length(files);

%% create directory to save results
if ~isdir(['results/' trackpars.title '/'])
    mkdir(['results/' trackpars.title '/']);
end







% clear
% clc
% 
% img_path = dir('shaking/*.jpg');
% gt = importdata(['shaking/' 'groundtruth.txt']);
% x = gt(:,1:2:end);
% y = gt(:,2:2:end);
% gt = [min(x,[],2), min(y,[],2), max(x,[],2) - min(x,[],2), max(y,[],2) - min(y,[],2)];
% 
% for i=1:size(img_path,1)
%     disp(i)
%     I = imread(['shaking/',img_path(i).name]);
%     figure    
%     imshow(I)
%     hold on
%     rectangle('Position', gt(i,:), 'linewidth', 3, 'EdgeColor', 'r');
%     imwrite(frame2im(getframe(gcf)),sprintf('../GT/shaking/%s.jpg', i));
%     close all
% end