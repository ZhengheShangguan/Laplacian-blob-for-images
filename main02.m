%% CS 543 Assignment 02
%% Zhenghe's Code
%% House Clean work
close all
clear all
clc

%% Initial work: Load the Image and Convert it into gray & double format
% name of the input file
imname = 'butterfly.jpg';
% read in the image
fullim = imread(imname);
% convert images to grayscale and double
fullim = im2double(rgb2gray(fullim));
% normalize the value to [0,1]
fullim = fullim/max(max(fullim));
% get the size of the image
[length, width] = size(fullim);

% Choose Laplacian filter method: 
% method 1: Increasing the kernel size by a factor k;
% method 2: Downsampling the image by a factor 1/k;
% method 3: Implement the difference-of-Gaussian pyramid.
method = 'method2';


%% Choose a method: Build the Laplacian scale space & Set Laplacian of Gaussian Filter
switch method
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'method1'
        %% Method 1: Increasing the kernel size by a factor k
        % Generate initial parameter values for Laplacian of Gaussian filter
        sigma = 1.5;
        Hsize = 5;
        k = 1.2;
        level = 15;
        thre = 1e-3;

        % initialize the scale space
        scale_space_method1 = zeros(length,width,level);
        
        % Nonmaximum suppression matrix, and its first and last layer will not be used
        image_max = zeros(length,width,level); 

        % tic for running time
        tic

        % For-loop for scale space
        for i = 1:1:level
            % generate a LoG filter and filter the image
            lap_gau = sigma^2 * fspecial('log', Hsize, sigma);
            scale_space_method1(:,:,i) = conv2(fullim, lap_gau, 'same').^2;

            % Save the corresponding radius for this layer
            radius(i) = sqrt(2)*sigma;

            % First, local nonmaximum suppression
            Lsize = ceil(Hsize*1.6);
            image_tmp(:,:,i) = ordfilt2(scale_space_method1(:,:,i),(Lsize)^2,ones(Lsize));
            % fun = @(x) max(max(x));
            % image_tmp(:,:,i) = nlfilter(scale_space_method1(:,:,i),[Lsize Lsize],fun);

            % Second, every single layer comparing to its neibors to get maximum for nonmax suppression
            if i>2
                tmp = cat(3,image_tmp(:,:,i-2),image_tmp(:,:,i-1),image_tmp(:,:,i));
                image_max(:,:,i-1) = max(tmp,[],3);                 
            end
            
            % Third, apply the threshold to original image_nonmax matrix for nonmax suppression
            scale_space_method1(:,:,i) = scale_space_method1(:,:,i).*(scale_space_method1(:,:,i)>thre);
            
            % Forth, update the size scale, and get the closest lower bound odd as size
            sigma = sigma * k;
            tmp = ceil(3*sigma);
            Hsize = tmp-mod(tmp+1,2);

            % toc for running time
            toc
        end
        
        % For-loop for neighbor's nonmaximum suppression
        MaxMap = cell(level,1);
        for i = 1:level
            % suppression by comparing with its two neighbors
            if (i>1)&&(i<level)
                Max_ind_matrix = (image_max(:,:,i)==scale_space_method1(:,:,i));
            elseif (i==1)
                Max_ind_matrix = (image_max(:,:,2)==scale_space_method1(:,:,i));
            else
                Max_ind_matrix = (image_max(:,:,level-1)==scale_space_method1(:,:,i));
            end
            MaxMap{i} = Max_ind_matrix;
            [row{i},col{i}] = find(Max_ind_matrix);
            num = sum(sum(Max_ind_matrix));
            rad{i} = radius(i).*ones(num,1);
        end
        row = vertcat(row{:});
        col = vertcat(col{:});
        rad = vertcat(rad{:});

        % plot
        show_all_circles(fullim, col, row, rad)


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'method2'
        %% Method 2: Downsampling the image by a factor 1/k
        % Generate parameter values for Laplacian of Gaussian filter
        sigma = 1.5;
        Hsize = 5;
        lap_gau = fspecial('log', Hsize, sigma);

        % Save the corresponding radius
        radius = sqrt(2)*sigma;

        % Parameters for scaling the image
        k = 1.2;
        level = 15;
        thre = 1e-4;

        % creates a cell array with level "slots"
        scale_space_method2 = cell(level,1); 

        % Nonmaximum suppression matrix, and its first and last layer will not be used
        image_max = zeros(length,width,level);        

        % tic for running time
        tic

        % For-loop for creating scale-space and local single layer nonmaximum suppression
        for i = 1:level
            
            % resize the image
            image_sc = imresize(fullim,1/(k^(i-1)));
            scale_space_method2{i} = conv2(image_sc, lap_gau, 'same').^2;
            
            % First, get the filtered matrix of the same size by resizing
            image_nonmax(:,:,i) = imresize(scale_space_method2{i}, [length, width], 'bicubic');
            
            % Second, local nonmaximum suppression in (Lsize x Lsize) local matrix
            Lsize = ceil(Hsize*k^(i-1));
            image_tmp(:,:,i) = ordfilt2(image_nonmax(:,:,i),Lsize^2,ones(Lsize));
            % fun = @(x) max(max(x));
            % mx = nlfilter(scale_space_method2(:,:,i),[Lsize,Lsize],fun);
     
            % Third, every single layer comparing to its neibors to get maximum for nonmax suppression
            if i>2
                tmp = cat(3,image_tmp(:,:,i-2),image_tmp(:,:,i-1),image_tmp(:,:,i));
                image_max(:,:,i-1) = max(tmp,[],3);                 
            end
            
            % Fourth, apply the threshold to original image_nonmax matrix for nonmax suppression
            image_nonmax(:,:,i) = image_nonmax(:,:,i).*(image_nonmax(:,:,i)>thre);

            % toc for running time
            toc
        end
        
        % For-loop for neighbor's nonmaximum suppression
        MaxMap = cell(level,1);
        for i = 1:level
            % suppression by comparing with its two neighbors
            if (i>1)&&(i<level)
                Max_ind_matrix = (image_max(:,:,i)==image_nonmax(:,:,i));
            elseif (i==1)
                Max_ind_matrix = (image_max(:,:,2)==image_nonmax(:,:,i));
            else
                Max_ind_matrix = (image_max(:,:,level-1)==image_nonmax(:,:,i));
            end
            MaxMap{i} = Max_ind_matrix;
            [row{i},col{i}] = find(Max_ind_matrix);
            num = sum(sum(Max_ind_matrix));
            rad{i} = (k^(i-1))*radius*ones(num,1);
        end
        
        % Combine the corresponding row, col, rad vectors
        row = vertcat(row{:});
        col = vertcat(col{:});
        rad = vertcat(rad{:});

        % plot
        show_all_circles(fullim, col, row, rad)
        
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'method3'
%% Bonus Part I: 
%% Implement the difference-of-Gaussian pyramid
        % Parameter values for scale space
        interval_num = 3;
        stack_num = interval_num + 3;
        k_Sigma = 2^(1/interval_num);
        SigmaOrigin = sqrt(2); % original sigma for the very first level
        Hsize = 21; % Gaussian filter size
        level = 3;
        
        % Scale Space Initialization
        Octave_Image = {level, stack_num};
        Octave_ImageDiff = {level, stack_num-1};
        radius = {level, stack_num};
        
        % tic for running time
        tic
        
        % For-loop for efficient SIFT detection with Gaussian Difference
        scale_space_method3 = cell(level,1);
        for Octave = 1:level
            % SigmaOrigin_level is the original sigma value for this level
            SigmaOrigin_level = SigmaOrigin * 2^(Octave-1);
            scale_space_method3{Octave} = imresize(fullim, 2^(-(Octave-1)));
            
            % do the convolution for all images in the same Octave level
            for i = 1:stack_num
                % set up the gaussian filter kernel
                Sigma_scale = SigmaOrigin_level * k_Sigma^(i-2);
                gau_filter = fspecial('gaussian', Hsize, Sigma_scale);
                % do the convolution
                Octave_Image{Octave, i} = conv2(scale_space_method3{Octave},gau_filter,'same').^2;
                        
                % save the radius
                radius{Octave,i} = sqrt(2)*Sigma_scale;
            end
            
            % toc for running time
            toc    

        end
        
        % Save the Gaussian difference as the result
        for Octave = 1:level
            for i = 1:stack_num-1
                Octave_ImageDiff{Octave,i} = Octave_Image{Octave, i+1} - Octave_Image{Octave, i};
            end
        end
        
        % Nonmaximum suppression for both neighbors of the pixel and the level
        MaxMap = cell(level,stack_num-3);
        % Initialize the row, col and rad matrix for saving the maximum pixel information
        row = {level, stack_num-2};
        col = {level, stack_num-2};
        rad = {level, stack_num-2};
        
        % for-loop for nonmax suppression
        for Octave = 1:level
            for i = 2:stack_num-2
                % initialize the maximum index matrix, mind the dimensions
                local_matrix_tmp = zeros(size(Octave_Image{Octave, i},1)-2, size(Octave_Image{Octave, i},2)-2, 3^3);
                ind = 0; % for counting
                
                % this for-loop for layers' nonmax suppression
                for j = i-1:i+1
                    % this for-loop for local pixel's nonmax suppression
                    for k = 1:3*3
                        % find out the subscripts of each pixel neighbors
                        [m,n] = ind2sub([3,3], k);
                        % save this local image
                        local_matrix_tmp(:,:,ind*9+k) = Octave_ImageDiff{Octave, j}(m:end-(3-m), n:end-(3-n));
                    end
                    ind = ind + 1;
                end
                
                % Nonmaximum suppression among the third dimension, we only
                % need the subscript matrix
                [~,Max_ind_matrix] = max(local_matrix_tmp, [], 3);
                % Only keep this layer's local maximum, whose index will 
                % always be 14 = 9 + 5(the middle layer's center pixel)
                Max_ind_matrix = ((Max_ind_matrix==14));
                % resize the matrix back to the original image size
                Max_ind_matrix = imresize(Max_ind_matrix,[length,width],'nearest');
                

                % find the index
                [row{Octave,i},col{Octave,i}] = find(Max_ind_matrix);
                num = sum(sum(Max_ind_matrix));
                rad{Octave,i} = 2^(Octave-1)*radius{Octave, i}*ones(num,1);
                MaxMap{Octave, i-1} = Max_ind_matrix;
            end
        end
        
        % Combine the corresponding row, col, rad vectors
        row = vertcat(row{:,2:stack_num-2});
        col = vertcat(col{:,2:stack_num-2});
        rad = vertcat(rad{:,2:stack_num-2});

        % plot
        show_all_circles(fullim, col, row, rad)
        
        
end



%% Bonus Part II: 
%% Eliminating edge response
if (method == 'method2')
    % threshold for eliminating edge response (according to the reference paper)
    thre_edge_gamma = 10;
    thre_edge = (thre_edge_gamma+1)^2 / thre_edge_gamma; 
    MaxMap_final = MaxMap; 

    % PS: for Gaussian, the MaxMap is cell(level, stack_num-3), thus here I
    % only implement this for method 2, it's very easy to change it into
    % the format for method 2 or 3
    for l = 1:level
        Map = MaxMap{l};
        % Scale-Space-Corner Index
        scale_space_corner_ind = find(Map); 

        if isempty(scale_space_corner_ind)
            continue;
        end
        for j = 1:size(scale_space_corner_ind,1)
            [Row,Col] = ind2sub([size(Map,1),size(Map,2)], scale_space_corner_ind(j));

            if (Row <= 1) || (Row >= size(Map,1)) || (Col <= 1) || (Col >= size(Map,2))
                MaxMap_final{l}(Row,Col) = 0; % discard out of matrix boundary
                continue;
            end

            % according to the paper, calculate the corresponding values to
            % reduce the edge response
            D_yy = image_nonmax(Row+1,Col,l) - 2*image_nonmax(Row,Col,l) + image_nonmax(Row-1,Col,l);
            D_xx = image_nonmax(Row,Col+1,l) - 2*image_nonmax(Row,Col,l) + image_nonmax(Row,Col-1,l);
            D_xy = image_nonmax(Row-1,Col+1,l) - image_nonmax(Row-1,Col-1,l) - image_nonmax(Row+1,Col+1,l) + image_nonmax(Row+1,Col-1,l);
            TrH = D_xx + D_yy;
            DetH = D_xx*D_yy - D_xy^2;
            if ((TrH^2 / DetH) >= thre_edge)
                MaxMap_final{l}(Row,Col) = 0; % discard unstable extrema
            end
        end
    end



    % For-loop for neighbor's nonmaximum suppression
    ROW1 = cell(level,1);
    COL1 = cell(level,1);
    RAD1 = cell(level,1);
    for i = 1:level
        [ROW1{i},COL1{i}] = find(MaxMap_final{i});
        num = sum(sum(MaxMap_final{i}));
        RAD1{i} = (k^(i-1))*radius*ones(num,1);
    end

    % Combine the corresponding row, col, rad vectors
    ROW1 = vertcat(ROW1{:});
    COL1 = vertcat(COL1{:});
    RAD1 = vertcat(RAD1{:});

    % plot
    figure(2);
    show_all_circles(fullim, COL1, ROW1, RAD1)
end

