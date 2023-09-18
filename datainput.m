
time = acquired.data.time;

stimulus_info_1 = acquired.stim(1,1).data;
stimulus_info_2 = acquired.stim(1,2).data;
stimulus_info_3 = acquired.stim(1,3).data;
stimulus_sequence_1 = zeros(length(time), 1);
stimulus_sequence_2 = zeros(length(time), 1);
stimulus_sequence_3 = zeros(length(time), 1);
for i = 1:size(stimulus_info_1, 1)
    start_time = stimulus_info_1(i, 1);
    duration = stimulus_info_1(i, 2);   
    end_time = start_time + duration;
    
 
    start_idx = find(time == start_time, 1);
    end_idx = find(time == end_time, 1);
    
    if isempty(start_idx) || isempty(end_idx)
        disp(['Warning: Start or end time not found in time series mapping']);
        continue;
    end
    
   
    stimulus_sequence_1(start_idx:end_idx) = 1;
end


for i = 1:size(stimulus_info_2, 1)
    start_time = stimulus_info_2(i, 1);
    duration = stimulus_info_2(i, 2);  
    end_time = start_time + duration;
    
   
    start_idx = find(time == start_time, 1);
    end_idx = find(time == end_time, 1);
    
    if isempty(start_idx) || isempty(end_idx)
        disp(['Warning: Start or end time not found in time series mapping']);
        continue;
    end
    
   
    stimulus_sequence_2(start_idx:end_idx) = 1;
end


for i = 1:size(stimulus_info_3, 1)
    start_time = stimulus_info_3(i, 1); 
    duration = stimulus_info_3(i, 2);  
    end_time = start_time + duration;
    
  
    start_idx = find(time == start_time, 1);
    end_idx = find(time == end_time, 1);
    
    if isempty(start_idx) || isempty(end_idx)
        disp(['Warning: Start or end time not found in time series mapping']);
        continue;
    end
    
    
    stimulus_sequence_3(start_idx:end_idx) = 1;
end

hbO_data = outResults.dc.dataTimeSeries(:, 1:61);


max_duration = min(30, max(time));  

tau_p = 2;
tau_d = 5;
A = 6;
t = time(time <= max_duration);  

hrf = (t.^tau_p).*exp(-t)/factorial(tau_p) - (t.^(tau_p + tau_d).*exp(-t))/(A*factorial(tau_p + tau_d));

figure;
plot(t, hrf, 'b-', 'LineWidth', 2);
xlabel('时间 (秒)');
ylabel('hrf');
title('Hemodynamic Response Function (HRF)');
grid on;




convolution_result_1 = conv(hrf, stimulus_sequence_1);
convolution_result_2 = conv(hrf, stimulus_sequence_2);
convolution_result_3 = conv(hrf, stimulus_sequence_3);
num_channels = 61;
betas_hbO_S1 = zeros(num_channels, 2);
betas_hbO_S2 = zeros(num_channels, 2);
betas_hbO_S3 = zeros(num_channels, 2);

for i = 1:num_channels
    X_HbO2_1 = [convolution_result_1, ones(length(convolution_result_1), 1)];
    X_HbO2_2 = [convolution_result_2, ones(length(convolution_result_2), 1)];
    X_HbO2_3 = [convolution_result_3, ones(length(convolution_result_3), 1)];
    X_HbO2_1 = X_HbO2_1(1:length(hbO_data), :);
    X_HbO2_2 = X_HbO2_2(1:length(hbO_data), :);
    X_HbO2_3 = X_HbO2_3(1:length(hbO_data), :);
    betas_hbO_S1(i, :) = X_HbO2_1 \ hbO_data(:, i); 
    betas_hbO_S2(i, :) = X_HbO2_2 \ hbO_data(:, i);
    betas_hbO_S3(i, :) = X_HbO2_3 \ hbO_data(:, i);
end



% Intercepts elements with the same number of rows as hbO_data
convolution_result_subset = convolution_result_3(1:size(hbO_data, 1), :);

% Perform element-level division
slope_1 = hbO_data ./ convolution_result_subset;

% slope_1 now contains the result of dividing the elements of the two matrices, with each column corresponding to one element of the result of the division

% Create an array of cells to store the data in slope_1.
slope_data_storage = cell(1, 5);

for i = 1:size(stimulus_info_3, 1)
    start_time = stimulus_info_3(i, 1); 
    duration = stimulus_info_3(i, 2);   
    end_time = start_time + duration;
    
    % Find the index corresponding to the start time and end time
    start_idx = find(time == start_time, 1);
    end_idx = find(time == end_time, 1);
    
    if isempty(start_idx) || isempty(end_idx)
        disp(['Warning: Start or end time not found in time series mapping']);
        continue;
    end
    
    % Extract the row intervals in slope_1 and store them in the cell array
    slope_data_storage{i} = slope_1(start_idx+1:end_idx, :);
end

% Now, slope_data_storage contains 5 different row intervals of data

% Create an array of cells to store the results of the t-tests
ttest_results = cell(1, 5);

for i = 1:5
    % Get current matrix from slope_data_storage
    current_matrix = slope_data_storage{i};
    
    % Perform a t-test assuming a mean of 0
    [~, p_values] = ttest(current_matrix, 0, 'Dim', 1);  % 在每列上执行 t-检验
    
    % Store p values in ttest_results
    ttest_results{i} = p_values;
end

% ttest_results contains the t-test results for each matrix

% Create a matrix of size 5x61, initialized to NaN
heatmap_data = NaN(5, 61);

%Populate the t-test results into the matrix
for i = 1:5
    ttest_result = ttest_results{i};
    heatmap_data(i, 1:length(ttest_result)) = ttest_result;
end

figure;
imagesc(heatmap_data);

colorbar;

title('T-Test Results Heatmap');
xlabel('T-Test Index');
ylabel('Array Index');



beta_values1 = betas_hbO_S1(:, 1);
beta_values2 = betas_hbO_S2(:, 1);
beta_values3 = betas_hbO_S3(:, 1);

zero_indices = find(plotOptions.blockActive{1, 1}(:, 3) == 0);


for i = 1:numel(zero_indices)
    if zero_indices(i) > 61
        zero_indices(i) = zero_indices(i) - 61;
    end
end
beta_values1(zero_indices) = [];
beta_values2(zero_indices) = [];
beta_values3(zero_indices) = [];


