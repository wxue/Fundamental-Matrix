% Matlab R2013a

clc; clear all; close all;

%% Settings
Threshold_R = 1/1000;               % threshold for corner response R
sigma = 1;
k = 0.12;
Threshold_C = 0.87;                 % threshold for correlation
RANSAC_8point_times = 10;
Threshold_RANSAC_eF = 0.1;
Threshold_RANSAC_selectpoints = 0.1;
Threshold_RANSAC_F_final = 0.1;

% 
%% Load both images & detect corners to extract patches
%
for pic = 1:2
    if pic == 1
        IMG = imread('uttower_left.jpg');       % read source JPEG image into 3D array
    else 
        IMG = imread('uttower_right.jpg');
    end
    IMG_hsv = rgb2hsv(IMG);                     % convert image from RGB to HSV

    IMG_hsv_value = IMG_hsv(:,:,3);             % gray scale value
    [Ix, Iy] = gradient(IMG_hsv_value);         % Shifts

    Ix2 = Ix.^2;
    Iy2 = Iy.^2;
    Ixy = Ix.*Iy;

    % Window Function
    windowsize = 3*sigma;
    [Wx, Wy] = meshgrid(-windowsize : windowsize, -windowsize : windowsize);
    w = exp(-(Wx .^ 2 + Wy .^ 2) / (2 * sigma ^ 2));

    % Convolutions
    A = conv2(w, Ix2);
    B = conv2(w, Iy2);
    C = conv2(w, Ixy);

    [x,y] = size(IMG_hsv_value);                % x,y is the width and length of the image
    R = zeros(x, y);                            % Initialize Corner response R
    % Loop for each pixel
    for xRows = 1:x
        for yColumns = 1:y  
            M = [A(xRows, yColumns), C(xRows, yColumns); C(xRows, yColumns), B(xRows, yColumns)];   % get M at each pixel
            R(xRows,yColumns) = det(M) - k * (trace(M) ^ 2);                                        % compute R at each pixel
        end
    end

    IMG_result = IMG;

    count = 1;
    for xRows = 3:x-2
        for yColumns = 3:y-2
            Local = [R(xRows-1, yColumns), R(xRows+1, yColumns), R(xRows, yColumns-1), R(xRows, yColumns+1)];
            if ((R(xRows, yColumns) > Threshold_R) && (R(xRows, yColumns) > max(Local)))
            % For those corner response R larger than Threshold and also local maximum
                IMG_result(xRows, yColumns, :) = [0, 0, 255];       % Mark corner point by blue
                if pic == 1
                    IMG_Left_corners(count, 1) = xRows;
                    IMG_Left_corners(count, 2) = yColumns;
                else 
                    IMG_Right_corners(count, 1) = xRows;
                    IMG_Right_corners(count, 2) = yColumns;
                end
                for xPatche = -2:2
                    for yPatche = -2:2 
                        if pic == 1 
                            Patches_left(count, xPatche+3, yPatche+3) = IMG_hsv_value(xRows+xPatche, yColumns+yPatche);
                        else
                            Patches_right(count, xPatche+3, yPatche+3) = IMG_hsv_value(xRows+xPatche, yColumns+yPatche);
                        end
                    end
                end
                count = count + 1;
            end
        end
    end
    
    if pic == 1
        corner_count_left = count-1;                                % Corner deceted mount
    else
        conrner_count_right = count-1;
    end
    
    % show results
    figure('Name', 'Corner Detected');  
    imshow(IMG_result);
end

%
%% Computing correlation & selcet top
%
count = 1;
for xLeft = 1:corner_count_left
    for yRight = 1:conrner_count_right
        c_value = corr2(Patches_left(xLeft,:,:), Patches_right(yRight,:,:));
        if c_value > Threshold_C
        % Record those correlation larger than Threshold(top)
            correlation(count, 1) = xLeft;      % corner point in the left image(in raster scan order)
            correlation(count, 2) = yRight;     % corner point in the right image(in raster scan order)
            correlation(count, 3) = c_value;    % correlation value
            count = count + 1
        end
    end
end
correlation_pairs_count = count -1;

count = 1;
for pair = 1:correlation_pairs_count            % loop each pair
    max_flag = 1;
    for check = 1:correlation_pairs_count       % make sure one point only be paired once
        if (correlation(pair, 1) == correlation(check, 1))&&(correlation(pair, 3) < correlation(check, 3))
            max_flag = 0;
            break
        end
        if (correlation(pair, 2) == correlation(check, 2))&&(correlation(pair, 3) < correlation(check, 3))
            max_flag = 0;
            break
        end
    end
    if max_flag
        correlation_pairs(count, 1) = correlation(pair, 1);
        correlation_pairs(count, 2) = correlation(pair, 2);
        correlation_pairs(count, 3) = correlation(pair, 3);
        count = count + 1;
        correlation_pairs
    end
end
correlation_pairs_count = count - 1;

%
%% 8-point RANSAC
%
for RANSAC_8point = 1:RANSAC_8point_times
    Rand_points = randi([1 correlation_pairs_count],1,8);
    for count = 1:8
        x_o(count,:) = IMG_Left_corners(correlation_pairs(Rand_points(count), 1), :);
        x_p(count,:) = IMG_Right_corners(correlation_pairs(Rand_points(count), 2), :);
    end
    UV = [x_o(1,1)*x_p(1,1), x_o(1,1)*x_p(1,2), x_o(1,1), x_o(1,2)*x_p(1,1), x_o(1,2)*x_p(1,2), x_o(1,2), x_p(1,1), x_p(1,2);
          x_o(2,1)*x_p(2,1), x_o(2,1)*x_p(2,2), x_o(2,1), x_o(2,2)*x_p(2,1), x_o(2,2)*x_p(2,2), x_o(2,2), x_p(2,1), x_p(2,2);
          x_o(3,1)*x_p(3,1), x_o(3,1)*x_p(3,2), x_o(3,1), x_o(3,2)*x_p(3,1), x_o(3,2)*x_p(3,2), x_o(3,2), x_p(3,1), x_p(3,2);
          x_o(4,1)*x_p(4,1), x_o(4,1)*x_p(4,2), x_o(4,1), x_o(4,2)*x_p(4,1), x_o(4,2)*x_p(4,2), x_o(4,2), x_p(4,1), x_p(4,2);
          x_o(5,1)*x_p(5,1), x_o(5,1)*x_p(5,2), x_o(5,1), x_o(5,2)*x_p(5,1), x_o(5,2)*x_p(5,2), x_o(5,2), x_p(5,1), x_p(5,2);
          x_o(6,1)*x_p(6,1), x_o(6,1)*x_p(6,2), x_o(6,1), x_o(6,2)*x_p(6,1), x_o(6,2)*x_p(6,2), x_o(6,2), x_p(6,1), x_p(6,2);
          x_o(7,1)*x_p(7,1), x_o(7,1)*x_p(7,2), x_o(7,1), x_o(7,2)*x_p(7,1), x_o(7,2)*x_p(7,2), x_o(7,2), x_p(7,1), x_p(7,2);
          x_o(8,1)*x_p(8,1), x_o(8,1)*x_p(8,2), x_o(8,1), x_o(8,2)*x_p(8,1), x_o(8,2)*x_p(8,2), x_o(8,2), x_p(8,1), x_p(8,2)];
    F = UV\[-1,-1,-1,-1,-1,-1,-1,-1]';
    F_matrix{RANSAC_8point} = [F(1), F(2), F(3); F(4), F(5), F(6); F(7), F(8), 1];
    F_matrix_weight{RANSAC_8point} = 0;
    for count = 1:correlation_pairs_count
        Xt = [IMG_Left_corners(correlation_pairs(count, 1), 1), IMG_Left_corners(correlation_pairs(count, 1), 2), 1];    
        Xp = [IMG_Right_corners(correlation_pairs(count, 2), 1); IMG_Right_corners(correlation_pairs(count, 2), 2); 1];
        if abs(Xt*F_matrix{RANSAC_8point}*Xp) < Threshold_RANSAC_eF
            F_matrix_weight{RANSAC_8point} = F_matrix_weight{RANSAC_8point} + 1;
        end
    end
end
M_value = 0;
for count = 1:RANSAC_8point_times               % Find a F_Matrix which fixes most
    if F_matrix_weight{count} > M_value
        M_value = F_matrix_weight{count};
        M_num = count;
    end
end

selcet_pairs_count = 1;
for count = 1:correlation_pairs_count
    Xt = [IMG_Left_corners(correlation_pairs(count, 1), 1), IMG_Left_corners(correlation_pairs(count, 1), 2), 1];    
    Xp = [IMG_Right_corners(correlation_pairs(count, 2), 1); IMG_Right_corners(correlation_pairs(count, 2), 2); 1];
    if abs(Xt*F_matrix{M_num}*Xp) < Threshold_RANSAC_selectpoints
        selcet_pairs(selcet_pairs_count, :) = correlation_pairs(count, :);
        selcet_pairs_count = selcet_pairs_count + 1;
    end
end
selcet_pairs_count = selcet_pairs_count - 1;

%
%% 8-point RANSAC with select pairs
%
for RANSAC_8point = 1:RANSAC_8point_times
    Rand_points = randi([1 selcet_pairs_count],1,8);
    for count = 1:8
        x_o(count,:) = IMG_Left_corners(selcet_pairs(Rand_points(count), 1), :);
        x_p(count,:) = IMG_Right_corners(selcet_pairs(Rand_points(count), 2), :);
    end
    UV = [x_o(1,1)*x_p(1,1), x_o(1,1)*x_p(1,2), x_o(1,1), x_o(1,2)*x_p(1,1), x_o(1,2)*x_p(1,2), x_o(1,2), x_p(1,1), x_p(1,2);
          x_o(2,1)*x_p(2,1), x_o(2,1)*x_p(2,2), x_o(2,1), x_o(2,2)*x_p(2,1), x_o(2,2)*x_p(2,2), x_o(2,2), x_p(2,1), x_p(2,2);
          x_o(3,1)*x_p(3,1), x_o(3,1)*x_p(3,2), x_o(3,1), x_o(3,2)*x_p(3,1), x_o(3,2)*x_p(3,2), x_o(3,2), x_p(3,1), x_p(3,2);
          x_o(4,1)*x_p(4,1), x_o(4,1)*x_p(4,2), x_o(4,1), x_o(4,2)*x_p(4,1), x_o(4,2)*x_p(4,2), x_o(4,2), x_p(4,1), x_p(4,2);
          x_o(5,1)*x_p(5,1), x_o(5,1)*x_p(5,2), x_o(5,1), x_o(5,2)*x_p(5,1), x_o(5,2)*x_p(5,2), x_o(5,2), x_p(5,1), x_p(5,2);
          x_o(6,1)*x_p(6,1), x_o(6,1)*x_p(6,2), x_o(6,1), x_o(6,2)*x_p(6,1), x_o(6,2)*x_p(6,2), x_o(6,2), x_p(6,1), x_p(6,2);
          x_o(7,1)*x_p(7,1), x_o(7,1)*x_p(7,2), x_o(7,1), x_o(7,2)*x_p(7,1), x_o(7,2)*x_p(7,2), x_o(7,2), x_p(7,1), x_p(7,2);
          x_o(8,1)*x_p(8,1), x_o(8,1)*x_p(8,2), x_o(8,1), x_o(8,2)*x_p(8,1), x_o(8,2)*x_p(8,2), x_o(8,2), x_p(8,1), x_p(8,2)];
    F = UV\[-1,-1,-1,-1,-1,-1,-1,-1]';
    F_matrix{RANSAC_8point} = [F(1), F(2), F(3); F(4), F(5), F(6); F(7), F(8), 1];
    F_matrix_weight{RANSAC_8point} = 0;
    for count = 1:selcet_pairs_count
        Xt = [IMG_Left_corners(selcet_pairs(count, 1), 1), IMG_Left_corners(selcet_pairs(count, 1), 2), 1];    
        Xp = [IMG_Right_corners(selcet_pairs(count, 2), 1); IMG_Right_corners(selcet_pairs(count, 2), 2); 1];
        if abs(Xt*F_matrix{RANSAC_8point}*Xp) < Threshold_RANSAC_F_final
            F_matrix_weight{RANSAC_8point} = F_matrix_weight{RANSAC_8point} + 1;
        end
    end
end
M_value = 0;
for count = 1:RANSAC_8point_times               % Find the result F_Matrix which fixes most
    if F_matrix_weight{count} > M_value
        M_value = F_matrix_weight{count};
        M_num = count;
    end
end

F_matrix_final = F_matrix{M_num}