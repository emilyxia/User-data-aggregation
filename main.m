% script to analysis StallCatchers.com data
% Developed by Mohammad Haft-Javaherian (mh973@cornell.edu)
% Fall 2017

%% Reading webiste data
% D columns = [date,user_id,movie_id,answer,ground_truth,type_train]
D = readtable('2017-09-12_SC_annotations.csv');
D = [exceltime(datetime(table2cell(D(:,1)), 'InputFormat', ...
    'yyyy-MM-dd HH:mm:ss')), table2array(D(:, 2:end))];


%% Adjacency and Accuracy metrics
% Initialization  
numUser = max(D(:, 2));
numVessel = max(D(:, 3));
A = NaN([numVessel, numUser]);
user_accuracy_hist = cell([numVessel, numUser]);
train = false(numVessel, 1);
ground_truth = NaN(numVessel, 1);

% Clean data
D(D(:, 4)==3, :) = [];
[~, i] = sort(D(:, 1));
D = D(i,:);

% Generating look up tables:
%   response code: 0=flow 1=stall NA=not available,
%   A: adjacency with [vessels,users],  
%   train: vessel ID of training vessels,
%   ground_truth: ground truth of 
%   user_accuracy_hist: history of each user prior to the current answer.  
for i=1:size(D,1)
    if mod(i,100000)==0,fprintf('%d-',i),end
    if isnan(A(D(i, 3),D(i, 2)))
        A(D(i, 3),D(i, 2)) = D(i, 4);
        user_accuracy_hist{D(i, 3),D(i, 2)} = D(find(D(1:i-1,2)==D(i,2).*D(1:i-1,6)==1),:); 
        if D(i, 6) == 1
            train(D(i, 3)) = true;
        end
        if D(i, 5) ~= 2
            ground_truth(D(i, 3)) = D(i, 5);
        end
    end
end
fprintf('\n',i)
user_accuracy_a_backup = user_accuracy_hist;

%% userAccuracy based on all training

% userAccuracy column : [TP, TN, FP, FN, sensitivity, specificity, dPrime,
% Jaccard, Dice]
userAccuracy = NaN(numUser, 9);
train = find(train);
test = unique(D(and(D(:, 5)<2, D(:, 6)==0), 3));
pd = makedist('Normal', 0, 1);
train_data = false(size(D,1),1);
train_data(D(:, 6) == 1) = true;

for i=1:numUser
    temp = D(and(D(:, 2) == i, train_data), 4:5);
    n = size(temp, 1);
    userAccuracy(i, 1) = sum((temp(:,1)).*(temp(:,2))) / n; % TP
    userAccuracy(i, 2) = sum((1-temp(:,1)).*(1-temp(:,2))) / n; % TN
    userAccuracy(i, 3) = sum((temp(:,1)).*(1-temp(:,2))) / n; % FP
    userAccuracy(i, 4) = sum((1-temp(:,1)).*(temp(:,2))) / n; % FN
    userAccuracy(i, 5) = userAccuracy(i, 1) / (userAccuracy(i, 1) + ...
        userAccuracy(i, 4)); % sensitivity
    userAccuracy(i, 6) = userAccuracy(i, 2) / (userAccuracy(i, 2) + ...
        userAccuracy(i, 3)); % specificity
    userAccuracy(i, 7) = max([0, icdf(pd, userAccuracy(i, 5)) - ...
        icdf(pd, 1-userAccuracy(i, 6))]); % dPrime
    userAccuracy(i, 8) = userAccuracy(i, 1) / (userAccuracy(i, 1) + ...
        userAccuracy(i, 3)+userAccuracy(i, 4)); % Jaccard
    userAccuracy(i, 9) = 2*userAccuracy(i, 1) / (2*userAccuracy(i, 1) + ...
        userAccuracy(i, 3)+userAccuracy(i, 4)); % Dice
end

%% calculating results based on all the methods
test_data = train;
A_train = A(train,:);
figure
hold on
AUC =[];
for k=1:8
    switch k
        case 1
            % sensitivity
            W = sum(userAccuracy(:, 5)'.*A, 2, 'omitnan') ./ ...
                sum(userAccuracy(:, 5)'.*(~isnan(A)), 2, 'omitnan');
        case 2
            % specificity
            W = 1-sum(userAccuracy(:, 6)'.*(1-A), 2, 'omitnan') ./ ...
                sum(userAccuracy(:, 6)'.*(~isnan(A)), 2, 'omitnan');
        case 3
            % dPrime
            W = sum(userAccuracy(:, 7)'.*(A), 2, 'omitnan') ./ ...
                sum(userAccuracy(:, 7)'.*(~isnan(A)), 2, 'omitnan');
        case 4
            % Jaccard
            W = sum(userAccuracy(:, 8)'.*(A), 2, 'omitnan') ./ ...
                sum(userAccuracy(:, 8)'.*(~isnan(A)), 2, 'omitnan');
        case 5
            % Dice
            W = sum(userAccuracy(:, 9)'.*(A), 2, 'omitnan') ./ ...
                sum(userAccuracy(:, 9)'.*(~isnan(A)), 2, 'omitnan');
        case 6
            % sens+spec
            W_P = sum(userAccuracy(:, 5)'.*A, 2, 'omitnan') ./ ...
                sum(userAccuracy(:, 5)'.*(~isnan(A)), 2, 'omitnan');
            W_N = sum(userAccuracy(:, 6)'.*(1-A), 2, 'omitnan') ./ ...
                sum(userAccuracy(:, 6)'.*(~isnan(A)), 2, 'omitnan');
            W = W_P ./(W_P+W_N);
        case 7
            % STAPLE1
            [~, p, q] = STAPLE_Silvia(A_train);
            [W1] = STAPLE_FWD(A(test_data,:), p, q);
            W(test_data) = W1;
        case 8
            % STAPLE2
            [W2, p, q] = STAPLE_Silvia(A(test_data,:));
            W(test_data) = W2;
    end

    % ROC and AUC generation and plots
    sensitivity = [];
    specificity = [];
    Jaccard = 0;
    Dice = 0;
    for W_limit=linspace(min(W(:)),max(W(:)),20)
        W0 = W >= W_limit;
        n = sum(~isnan(W0(test_data)));
        TP = sum((W0(test_data)) .* (ground_truth(test_data)), 'omitnan') /n;
        TN = sum((1-W0(test_data)) .* (1-ground_truth(test_data)), 'omitnan') /n;
        FP = sum((W0(test_data)) .* (1-ground_truth(test_data)), 'omitnan') /n;
        FN = sum((1-W0(test_data)) .* (ground_truth(test_data)), 'omitnan') /n;
        sensitivity = [sensitivity; TP / (TP + FN)];
        specificity = [specificity; TN / (TN + FP)];
        Jaccard = max([Jaccard, TP /(TP + FN + FN)]); 
        Dice = max([Dice, 2*TP /(2*TP + FN + FN)]);
    end
    AUC = [AUC, sum((-specificity(1:end-1)+specificity(2:end)).* ...
        (sensitivity(1:end-1)+sensitivity(2:end))/2)];
    plot(1-specificity,sensitivity)
end
line([0 1],[0 1])
legend({'Sensitivity', 'Specificity', 'd''','Jaccard', 'Dice', ...
    'Sens+Spec','STAPLE 1','STAPLE 2','y=x'})
ylabel('True Positive Rate') 
xlabel('False Positive Rate')
ylim([0 1])
xlim([0 1])
hold off
shg

% AUC plot
figure
bar(AUC)
ylabel('AUC')
ylim([0 1])
ax = gca;
ax.XTick = [1 2 3 4 5 6 7 8]; 
ax.XTickLabels = {'Sensitivity', 'Specificity', 'd''','Jaccard', 'Dice', ...
    'Sens+Spec','STAPLE 1','STAPLE 2'};
ax.XTickLabelRotation = 45;
shg


%% Time based userAccuracy 
% generate hitory look up tables just for test data
% this is for cased we do not want to generate it for all the record in D
train = find(train);
test = unique(D(and(D(:, 5)<2, D(:, 6)==0), 3));
train_data = false(size(D,1),1);
train_data(D(:, 6) == 1) = true;

test_data = test;
for j = test_data'
        for k = find(~isnan(A(j,:)))
            i = find((D(:, 3)==j).*(D(:, 2)==k),1);
            user_accuracy_hist{D(i, 3),D(i, 2)} = D(find(D(1:i-1,2)==D(i,2).*D(1:i-1,6)==1),:); 
        end
end

%%
figure
hold on
legend_text = {};
AUC = [];
isJustTestHistory = 1; % change this if the total history is available and we want to run the whole D

for  i = 1:25
    numHistoryRecord = 40 * i; 
    if isJustTestHistory
        JI = nan(size(A));
        for j = test_data'
            for k = find(~isnan(A(j,:)))
                temp = user_accuracy_hist{j,k}(max(1,end-numHistoryRecord+1):end,4:5);
                n = size(temp, 1);
                TP = sum((temp(:,1)).*(temp(:,2))) / n; % TP
                TN = sum((1-temp(:,1)).*(1-temp(:,2))) / n; % TN 
                JI(j,k) = TP /(1-TN);
            end
        end
    else
        TP = cellfun(@(x) sum(x(max(1,end-numHistoryRecord+1):end,4) .* ...
            x(max(1,end-numHistoryRecord+1):end,5)) / ...
            numel(x(max(1,end-numHistoryRecord+1):end,4)),user_accuracy_a);
        TN = cellfun(@(x) sum((1-x(max(1,end-numHistoryRecord+1):end,4)) .* ...
            (1-x(max(1,end-numHistoryRecord+1):end,5))) / ...
            numel(x(max(1,end-numHistoryRecord+1):end,4)),user_accuracy_a);
        JI = TP ./(1-TN);
    end     
    
    W = mean(JI.*A, 2, 'omitnan'); 
    
    % generate and plot ROC and AUC
    sensitivity = [];
    specificity = [];
    Jaccard = 0;
    Dice = 0;
    for W_limit=linspace(min(W(:)),max(W(:)),20)
        W0 = W >= W_limit;
        n = sum(~isnan(W0(test_data)));
        TP = sum((W0(test_data)) .* (ground_truth(test_data)), 'omitnan') /n;
        TN = sum((1-W0(test_data)) .* (1-ground_truth(test_data)), 'omitnan') /n;
        FP = sum((W0(test_data)) .* (1-ground_truth(test_data)), 'omitnan') /n;
        FN = sum((1-W0(test_data)) .* (ground_truth(test_data)), 'omitnan') /n;
        sensitivity = [sensitivity; TP / (TP + FN)];
        specificity = [specificity; TN / (TN + FP)];
        Jaccard = max([Jaccard, TP /(TP + FN + FN)]); 
        Dice = max([Dice, 2*TP /(2*TP + FN + FN)]);
    end
    AUC = [AUC, sum((-specificity(1:end-1)+specificity(2:end)).* ...
        (sensitivity(1:end-1)+sensitivity(2:end))/2)];
    plot(1-specificity,sensitivity)
    legend_text{i} = strcat(num2str(numHistoryRecord));
end
legend_text{i+1} = 'y=x';
line([0 1],[0 1])
legend(legend_text)
ylabel('True Positive Rate') 
xlabel('False Positive Rate')
ylim([0 1])
xlim([0 1])
hold off
shg
% plot AUC
figure
bar(AUC)
ylabel('AUC')
xlabel('Number of record considered for accuracy measurment')
ylim([floor(min(AUC*100)) ceil(max(AUC*100))]/100)

ax = gca;
legend_text = legend_text(1:end-1);
ax.XTick = 1:length(legend_text); 
ax.XTickLabels = legend_text;
ax.XTickLabelRotation = 45;
xlim([0 length(legend_text)+0.5])
shg


