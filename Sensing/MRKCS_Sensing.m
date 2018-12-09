function [R, G, results] =MRKCS_Sensing(imgSize, m , trial)
% Function to generate diagonal multiple resolution Gaussian sensing matrix
% with given ratio
%	[R, G, results] = MRKCS_Sensing(imgSize, subrate, num_res)
%	- Input:
%		+ imgSize: size of input image
%       + m: meaurement allocation vector 
%		+ trial: no. of trial for different sensing matrix
% 	- Output
%		+ R, G: KCS sensing matrix at imgSize resolution
%		+ results: contains all sub/low resolution

level = max(size(m)); 

% sensing matrix
R_all = cell(1); G_all = cell(1);
for i = 1:1:level-1
    if(i==1);
        [R_all{i}, G_all{i}] = KCS_SensingMtx(imgSize/(2^(level-i)), m(i), trial);
    end;
    [R_all{i+1}, G_all{i+1}] = KCS_SensingMtx(imgSize/(2^(level-i)), m(i+1), trial);
end;
Rtmp = R_all{1};
Gtmp = G_all{1};
for i =1:1:level-1
    [Rtmp, Gtmp]     = blk_KCS(Rtmp, R_all{i+1}, Gtmp, G_all{i+1}, imgSize/(2^(level-i-1)));
end;

R = Rtmp;
G = Gtmp;
% update additional info
results.R_all = R_all;
results.G_all = G_all;

results.subrate_real = size(R,1).^2/(imgSize.^2);

function [R, G] = KCS_SensingMtx(imgSize, m, trial)
% function to generate sensing matrix
patch = ['SensingMatrixSepTV' num2str(imgSize)];
if ~exist(patch, 'dir');
    display(['Generate new folder for sensing matrix']);
    mkdir(patch);
end;

projectionMatrixFile 	= [patch '\SenMtr_R_trial' num2str(trial) '.mat'];
R                       = SPL_GenerateProjection(imgSize, m, projectionMatrixFile);

projectionMatrixFile    = [patch '\SenMtr_G_trial' num2str(trial) '.mat'];
G                       = SPL_GenerateProjection(imgSize, m, projectionMatrixFile);
G                       = transpose(G);

function Phi = SPL_GenerateProjection(N, M, filename)

% N = block_size * block_size;
% M = round(sqrt(subrate) * N);

if ((nargin == 3) && exist(filename, 'file'))
    load(filename);
else
  Phi = orth(randn(N, N))';
end

if ((nargin == 3) && (~exist(filename, 'file')))
    display('Generate new sensing matrix')
    save(filename, 'Phi');
end

Phi = Phi(1:M, :);


function [R, G] = blk_KCS(R1, R2, G1, G2, img_size)
m1                    = size(R1, 1);
m2 					  = size(R2, 1);
R                     = zeros(m1 + m2, img_size);
R(1:m1, 1:img_size/2)        = R1;
R(m1+1:end,img_size/2+1:end) = R2;
m1                    = size(G1, 2);
m2 					  = size(G2, 2);
G                     = zeros(img_size, m1 + m2);
G(1:img_size/2, 1:m1)         = G1;
G(img_size/2+1:end, m1+1:end) = G2;