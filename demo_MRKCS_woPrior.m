
close all;  clear all; % clc;
path(path,genpath(pwd));

%% General setting
par.imgSize     = 512;
par.sparsity    = [0.1 ];
testIm          = [1 ];
par.testIm      = testIm; 
par.recMode     = 'DTVNL1'; 
par.filterMode  = 'BM3D';   % 'No', 'NLM', 'BM3D', 'WNNM'
par.sigma       = 0.02;    
par.showPSNR    = true;  % true/false
par.nbrLoop     = 1;

par.maxIter     = 25; 
par.qid_index   = [1 2 3  5:19];
par.level       = 3; 
par.tap         = 2; 
par.tau         = 10; 
par.aSeed       = 4; 
custorm_note    = ['_160127_tap' num2str(par.tap) '_ls' num2str(par.level)...
                   '_aSeed' num2str(par.aSeed)]; 

par.curLevel    = 4;   
recSize         = par.imgSize/(2^(par.level - par.curLevel +1));       % the size of recovered image            

for imgId = 1:1:length(testIm)
    [imgOrg, imgName] = testImage(par.imgSize, testIm(imgId)); 
    par.imgOrg = imgOrg; par.imgName = imgName;
    saveFolderText   = ['Result_text' num2str(par.imgSize) '\' ];    if ~exist(saveFolderText, 'dir'); mkdir(saveFolderText);   end;    
    fileNameNave  	 = [saveFolderText imgName '_' par.recMode '_' par.filterMode custorm_note ];
    write_info([fileNameNave  '.txt'], [par.imgName ' size' num2str(par.imgSize)]);
    write_info([fileNameNave  '.txt'], ['               size sigma   iter  sub    time  mse     psnr    ssim   mssim  vsnr   vif    vifp   uqi    ifc    nqm   snr    HFEN    srsim  rfsim  gmsd   fsim   psnrR  psnrG  psnrB psnrAvg psnrSt ssimR  ssimG  ssimB  ssimAv ssim   psnr1 psnr2 psnr3 ssim1 ssim2 ssim3']);

    resutls_all      = cell(1); 
    
    for sub = 1:1:length(par.sparsity)
        subrate = par.sparsity(sub);   
        
        qid_inter = cell(1);
        recImgAll = cell(1);
        for trial = 1:1:par.nbrLoop
            display(['Recover ' par.imgName ' using ' par.recMode  ' ' par.filterMode ', subrate' num2str(subrate) ',trial:' num2str(trial) '/' num2str(par.nbrLoop)]);
            
            % -------------- Sensing matrix ------------------------
            % Separabel wavelet transform 
            [W, W_all]= wavelet_matrix(par.imgSize, par.tap, par.level);
            
            % allocate measurement
            [m, rSubband] = measAlloc(subrate, par.imgSize, par.level, par.aSeed);
            par.measAllocate = m;  par.subrateAllocate = rSubband;          
                        
            % making sensing matrix			
            [R, G, results] = MRKCS_Sensing(par.imgSize, m , trial);
            
            % CS measuremetn 
            Y        = R*W*imgOrg*W'*G;
            
            % --------------- Recover the image -----------------------
            % Extract corresponding measurement and sensing matrix
            mq      = sum(m(1:par.curLevel));             
            Yq      = Y(1:mq, 1:mq);
            % extract sensing matrix
            Rq      = [];   Gq = [];            
            for i = 1:1:par.curLevel
                Rq = blkdiag(Rq, results.R_all{i}); 
                Gq = blkdiag(Gq, results.G_all{i});
            end;                                     
            % extract wavelet
            [Wq, ~] = wavelet_matrix(recSize, par.tap, par.curLevel - 1);
            normal = sqrt(2)^(par.level - par.curLevel + 1); 
            
            % recovery 
            par.imgOrg  = imresize(imgOrg, recSize/size(imgOrg, 1)); 
            tic;
            if size(Yq, 1) == size(Rq, 2)  % subrate = 1
                Wq  = Wq/normal; 
                recImg = (Rq*Wq)' *Yq * (Wq'*Gq)';
            else 
                Wq          = Wq * normal;
                opts        = setParams(par.recMode, par);                
                [recImg, ~] = DecWTVNLR(Rq*Wq, Wq'*Gq, Yq, opts, par.imgOrg, par.recMode);   
            end;     
            
            recTime(trial) = toc;      
            recImgAll{trial} = recImg; 
            
            % save each trial
            qid_inter{trial} = qid_cal(recImg, par.imgOrg, par.qid_index); 
            %write_results([fileNameNave  '.txt'], par.imgOrg, subrate, 0, qid_inter{trial}, 0, trial, par.nbrLoop);
            
            %results.resInter  = resInter;
            results.qid_inter = qid_inter;		
            results.t_org     = recTime;
            results.recImgAll = recImg;
            
            % save the results
            patch3 = ['Results\' par.recMode '_' par.filterMode '_iter' num2str(par.maxIter) '_sub' num2str(subrate*100) ];  if ~exist(patch3, 'dir');  mkdir(patch3); end;
            patch31 = [patch3 '\' par.imgName '_' par.recMode custorm_note '_trial' num2str(trial) '.mat']; save(patch31, 'results', 'par');
            if trial == 1
            save_image_result(recImg, custorm_note, imgName, subrate, qid_inter{trial}.psnr, qid_inter{trial}.ssim, trial);
            end;
        end; % end trial
        % save average
        write_results([fileNameNave  '.txt'], imgOrg, subrate, 0, qid_inter, 0, 0, trial);                
        display(['========== Recovery PSNR:' num2str(qid_inter{trial}.psnr) '============']);
       
        Note = [num2str(par.imgSize) '_' num2str(size(par.imgOrg,2)) '_sub' num2str(subrate) ];
        save_image_result(recImg, [ Note ] , imgName, subrate, qid_inter{trial}.psnr, qid_inter{trial}.ssim, 0);
            
		% calculate PSNR of block no boundary				
		results_all{sub} = results;         
        
    end; % end sparsity    
    % save all file 
    save([patch3 '\Full_' par.imgName '_' par.recMode custorm_note '.mat' ], 'results_all', 'par');
end; % end test image
display('END SIMULATION!!!');
