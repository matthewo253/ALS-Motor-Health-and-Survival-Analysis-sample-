clear;
clc;
close all;

% loads the als_clinical csv file
data = readtable('als_clinical.csv');
if(isempty(data))
	fprintf("No data table");
	return;
end

% initalizes arrays 
allPatientID = nan(1, numel(data.PatientID));
allTimeMonths = nan(1, numel(data.TimeMonths));
allALSFRS_R = nan(1, numel(data.ALSFRS_R));
allFVC = nan(1, numel(data.FVC));
allStatus = nan(1, numel(data.Status));

% puts values from the csv file into the arrays if they fit the parameters
n = 1;
for i = 1: numel(data.PatientID)
	if(~isnan(data.ALSFRS_R(i)) && ~isnan(data.FVC(i)))
		if((data.ALSFRS_R(i) < 48 && data.ALSFRS_R(i) > 0) && (data.FVC(i) < 100 && data.FVC(i) > 0))
			allPatientID(n) = data.PatientID(i);
			allTimeMonths(n) = data.TimeMonths(i);
			allALSFRS_R(n) = data.ALSFRS_R(i);
			allFVC(n) = data.FVC(i);
			allStatus(n) = data.Status(i);

			n = n + 1;
		end
	end
end

% sets parameters
maxALSFRS_R = 48;
motorNeuron = 10000;
respiratoryProxy = 8000;

% initalizes more arrays
motorHealthProxy = nan(1, numel(allALSFRS_R));
motorDegeneratedProxy = nan(1, numel(allALSFRS_R));
respiratoryHealthProxy = nan(1, numel(allFVC));
respiratoryDegeneratedProxy = nan(1, numel(allFVC));


% gets all the Healthy neurons and degenerated neurons
for i = 1: numel(allALSFRS_R)
	motorHealthProxy(i) = (allALSFRS_R(i) / maxALSFRS_R) * motorNeuron;
	motorDegeneratedProxy(i) = motorNeuron - motorHealthProxy(i);
end

for i = 1: numel(allFVC)
	respiratoryHealthProxy(i) = (allFVC(i) / 100) * respiratoryProxy;
	respiratoryDegeneratedProxy(i) = respiratoryProxy - respiratoryHealthProxy(i);
end

% creates arrays for the time span of healthy neurons
count = 0;
writeidx = 1;
allHealthyTimeMonths = nan(1, numel(allTimeMonths));
allHealthyMotorVals = nan(1, numel(allTimeMonths));
allHealthyPatientsID = nan(1, max(allTimeMonths));
for j = 1: numel(allTimeMonths) - 1
	for k = j + 1: numel(allTimeMonths)
		if(allTimeMonths(j) == allTimeMonths(k))
			count = count + 1;
		end

	end
	if(count >= 2)
		for offset = 0: count
            idx = j + offset;
            if(idx <= numel(allTimeMonths) && writeidx <= numel(allHealthyMotorVals))
			    allHealthyMotorVals(writeidx) = motorHealthProxy(idx);
			    allHealthyPatientsID(writeidx) = allTimeMonths(idx);
			    writeidx = writeidx + 1;
            end
		end
	end
	count = 0;
end




% initalizes arrays and creates the length of survival for the neurons
count = 0;
allHealthyBaseLine = nan(1, max(allPatientID));
survivalTime = nan(1, max(allPatientID));
censorFlags = ones(1, max(allPatientID));

m = 1;

cleanIDs = [];
for i = 1:length(allPatientID)
    if(~isnan(allPatientID(i)))
        cleanIDs(end+1) = allPatientID(i);
    end
end

uniqueIDs = unique(cleanIDs);

for i = 1:numel(uniqueIDs)
    currentID = uniqueIDs(i);
    indexes = [];
    for j = 1: numel(allPatientID)
	    if(allPatientID(j) == currentID)
		    indexes(end + 1) = j;
	    end
    end
    times = allTimeMonths(indexes);


    minIndex = 0;
    maxIndex = 0;

    maxVal = 0;
    minVal = 10000;

    for k = 1: numel(times)
	    if(~isnan(times(k)))
		    if(times(k) > maxVal) 
			    maxVal = times(k);
			    maxIndex = indexes(k);
		    end
    
		    if(times(k) < minVal)
			    minVal = times(k);
			    minIndex = indexes(k);
		    end
	    end
    end
	    
    if ~isnan(allALSFRS_R(minIndex))
        allHealthyBaseLine(m) = allALSFRS_R(minIndex);
        survivalTime(m) = allTimeMonths(maxIndex);
        if allStatus(maxIndex) == 1
            censorFlags(m) = 0;
        end
        m = m + 1;
    end
end


% creates an array of all the healthy patients	

validPatientVal = [];
for i = 1: numel(allHealthyBaseLine)
    if(~isnan(allHealthyBaseLine(i)))
        validPatientVal(end + 1) = allHealthyBaseLine(i); 
    end
end 

allSortedHealthyPatients = sort(validPatientVal);

numValidPatientVals = numel(allSortedHealthyPatients);

if mod(numValidPatientVals, 2) == 1
    medianMotorHealth = allSortedHealthyPatients((numValidPatientVals + 1) / 2);
else
    mid1 = allSortedHealthyPatients(numValidPatientVals / 2);
    mid2 = allSortedHealthyPatients(numValidPatientVals / 2 + 1);
    medianMotorHealth = (mid1 + mid2) / 2;
end

motorHealthDeviation = validPatientVal - medianMotorHealth;



% Creates a Histogram for the BaseLine Motor Health

figure;
histogram(allHealthyBaseLine, 20); % adjust bin number as needed
xlabel('Baseline Motor Health Proxy');
ylabel('Number of Patients');
title('Histogram of Baseline Motor Health Proxy');

% Scatter plot that shows ALSFRS_R scores vs the FVC scores
figure;
scatter(allALSFRS_R, allFVC, 20, 'filled');
xlabel('ALSFRS-R Score');
ylabel('FVC (%)');
title('ALSFRS-R vs. FVC');
grid on;


% Histogram of the Median deviation of the different health
figure;
histogram(motorHealthDeviation, 20);
xlabel('Deviation from Median Health');
ylabel('Number of Patients');
title('Deviation from Median Motor Health at Baseline');



% Kaplan-Meier Survival Curve

validIdx     = ~isnan(allHealthyBaseLine) & ~isnan(survivalTime);
baselineVals = allHealthyBaseLine(validIdx);
survivalVals = survivalTime(validIdx);
censorVals   = (censorFlags(validIdx)==1);
medVal = median(baselineVals, 'omitnan');
isHigh = baselineVals >= medVal;
isLow  = baselineVals <  medVal;

%fprintf('KM plot: total=%d, high=%d, low=%d, events=%d\n', ...
%    numel(baselineVals), sum(isHigh), sum(isLow), sum(~censorVals));
%disp('Survival durations (unique):'); disp(unique(survivalVals));

if ~isempty(baselineVals)
    [Sh,Xh] = ecdf(survivalVals(isHigh), ...
                   'Censoring', censorVals(isHigh), ...
                   'Function','survivor');

    Sh_new = zeros(length(Sh)+1, 1);
    Xh_new = zeros(length(Xh)+1, 1);
    Sh_new(1) = 1;
    Xh_new(1) = 0;
    for i = 1:length(Sh)
        Sh_new(i+1) = Sh(i);
        Xh_new(i+1) = Xh(i);
    end
    Sh = Sh_new;
    Xh = Xh_new;

    [Sl,Xl] = ecdf(survivalVals(isLow), ...
                   'Censoring', censorVals(isLow), ...
                   'Function','survivor');

    Sl_new = zeros(length(Sl)+1, 1);
    Xl_new = zeros(length(Xl)+1, 1);
    Sl_new(1) = 1;
    Xl_new(1) = 0;
    for i = 1:length(Sl)
        Sl_new(i+1) = Sl(i);
        Xl_new(i+1) = Xl(i);
    end
    Sl = Sl_new;
    Xl = Xl_new;

    % Creates the plot figures for the Kaplan Meier Survival Curve
    figure; hold on;
    stairs(Xh, Sh, 'r', 'LineWidth', 2);
    stairs(Xl, Sl, 'b', 'LineWidth', 2);
    xlabel('Time (Months)');
    ylabel('Survival Probability');
    legend('High Baseline','Low Baseline','Location','best');
    title('Kaplan–Meier Survival by Baseline Motor Proxy');
    grid on;
end
