function [W, W_all]= wavelet_matrix(img_size, tap, ds_level)
% [W, W_all]= wavelet_matrix(img_size, tap, ds_level)
%   Function to general separate wavelet amtrix, currently only support
%   up-to level 3;
%

W = eye(img_size);
W_all = cell(1);
for i = 1:1:ds_level
    tmp = wavelet_matrix_level(img_size, tap, i);
    W   = tmp*W;
    W_all{i} = tmp;
end

end

function W = wavelet_matrix_level(img_size, tap, ds_level)
% function to give wavelet matrix for each level of decomposition
%   W = wavelet_matrix(img_size, ds_level)
%   the wavelet coefficient of level 2 is given as
%       W_coef = W2*(W1* F * W1')*W2'.
%   The higher level of decomposition, can be esily follow.
%   - Input:
%       + img_size: size of input image
%       + tap: the specific wavelet in Daubeche family
%       + ds_level: level of decomposition
%   - Output: wavelet matrix
%       + W wavelet matrix for corresponding level

% Note: currently supprot level 2
switch ds_level
    case 1
        W = dau_mtx2(tap, img_size);
        
    case 2
        W1 = dau_mtx2(tap, img_size/2);
        W = make_block_based_wavelet(W1, img_size/2);
        
    case 3
        W1 = dau_mtx2(tap, img_size/4);
        W2 = make_block_based_wavelet(W1, img_size/4);
        W  = make_block_based_wavelet(W2, img_size/2);
    case 4
        W1 = dau_mtx2(tap, img_size/8);
        W2 = make_block_based_wavelet(W1, img_size/8);
        W3 = make_block_based_wavelet(W2, img_size/4);
        W  = make_block_based_wavelet(W3, img_size/2);
end
end

function W = make_block_based_wavelet(W1, img_size)
I  = eye(img_size);
[W, ~] = blk_KCS(W1, I, W1, I, img_size*2);
end

%% function to generate dau4 matrix
% daub4 tap
function W = dau_mtx2(tap, img_size)
% testing
% img_size = 512;
% tap = 6;

% daubeches 4 taps
scale_num = dbwavf(['db' num2str(tap/2)]) * sqrt(2);

wave_num = zeros(1, tap);
for i = 1:1:tap
    wave_num(i) = scale_num(tap-i+1)*( (-1 )^(i-1) );
end


% separate form
Wl = zeros(img_size/2, img_size);
Wh = zeros(img_size/2, img_size);
for i = 1:1:(img_size - tap)/2 +1
    if(i == 255)
        display('i = 255');
    end;
    z_bef = zeros(1, (i-1)*2);
    z_aft = zeros(1, img_size - (i-1)*2 - tap);
    Wl(i, :)   = [z_bef scale_num z_aft];
    Wh(i, :)   = [z_bef wave_num z_aft];
    
end;
iter = 0;
for j = i+1 : img_size/2
    iter        = iter + 1;
    index       = tap - iter*2 +1;
    Wl(j, :)    = [scale_num(index:end) zeros(1, img_size - tap) scale_num(1:index -1)];
    Wh(j, :)    = [wave_num(index:end) zeros(1, img_size - tap) wave_num(1:index -1)];
end;
W  = [Wl ; Wh];

%% display results
if(nargin > 2)
    testIm          = [1 ];
    img_size = 512;
    title_name = ['Wavelet Dau' num2str(tap)];
    if(tap ==2)
        title_name = ['Wavelet Haar'];
    end;
    
    for mmm = 1:length(testIm)
        [image, img_name] = testImage(img_size, testIm(mmm));
        wavelet_coef    = W*image*W';
        figure();
        imshow(wavelet_coef, []); title(title_name);
    end;
    figure();
    imshow(wavelet_coef(1:img_size/2, 1:img_size/2), []);
    title([title_name ' LL']);
end;

end

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
end