%% Rens Meerhoff
% 28-03-2019
% Script to generate mockdata for a paper that is in preparation with the
% title: "Exploring a Link Between Commitment to English Language Teaching
% and Student Learning Outcomes" by Moodie and Meerhoff.
%
% NB: Script works as a standalone, but then the teacher data is also
% simulated. For the paper, we used actual teacher data.

clear all
close all

%%
% Create a mock data set that represents the student pre and post test
% scores.

%% User Parameters
% Random numbers are generated along a normal distribution.
% The possibility exists to only accept generated numbers between a pre-defined minimum and/or maximum value.
% To make sure this is feasible, the standard deveiation must be roughly 1/8 of the range between the min and max parameters (assuming the average is in the middle).
%
% The parameters to determine the number of students per teacher.
% Note that the std should be roughly 1/8 of the range between the
% max and min parameters.
minimumNumberOfStudents = 4; % The script will re-generate the random numbers if one teacher has less students than this.
% NB: for student numbers there is no strict upper limit
nAverage = 30;
nStd = 10;

% The parameters to determine the pretest scores.
% Note that the std should be roughly 1/8 of the range between the
% max and min parameters.
pretestParameters.min = 5;
pretestParameters.max = 27-2; % taking max improvement into account
pretestParameters.avg = 15;
pretestParameters.std = 4;

% The parameters to determine the defaultImprovement.
% Note that the std should be roughly 1/8 of the range between the
% max and min parameters.
defaultImprovementParameters.min = -1.5;
defaultImprovementParameters.max = 4;
defaultImprovementParameters.avg = 1.25;
defaultImprovementParameters.std = 0.25;

% Bias in absolute points, where A_top gets on average a bias improvement,
% A_mid gets no improvement/decrement and A_bot gets on average a bias
% decrement.
% NB: Bias MUST BE smaller than or equal to 1
bias = [0.1 0.5 1.0];

%% / User Parameters

%% Input Data
% For the paper, we simulated the students' data, but we used real data
% from the teachers. In this public version of the script, we omitted the
% real teacher data and replaced it with made up data.

nTeachers = 70;

% Teacher IDs as found in "Full lists of pre-post tests non-anonymous.xlsx"
tID = 1 : nTeachers;

% Teacher scores
% 1 = Affective Field
% 2 = Continuance Field
% 3 = Normative Field
% 4 = Affective Work
% 5 = Continuance Work
% 6 = Normative Work
% 7 = Affective (summed)
% 8 = Continuance (summed)
% 9 = Normative (summed)
tScores = createRandomTeacherScores(nTeachers);

%% / Input Data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Generate a string to save the data as
fileString = ['RandomGen_' datestr(now, 'HH-MM-dd-mmm-yyyy')];

%% Student results - create mock data
% Generate the number of students per teacher
% nStudents = generateStudentNumbers(minimumNumberOfStudents, nAverage, nStd, tID);
nStudents = round(genSemiRandom(minimumNumberOfStudents, nAverage, nStd, length(tID)));

% Allocate the teachers to a rank based on total A score
groupsA = getGroupScores(7, tScores);
groupsC = getGroupScores(8, tScores);
groupsN = getGroupScores(9, tScores);

% Plot the student numbers
studentNumberBoxplots(nStudents, fileString, groupsA)

for i = 1:length(bias)
    %     % The old method
    %     generateMockData(bias(i), tID, nStudents, A_top, A_mid, A_bot, tScores);
    % The new method: normally distributed randomization
    generateMockData_normallyDistributed(bias(i), tID, nStudents, groupsA, groupsC, groupsN, tScores, fileString, pretestParameters, defaultImprovementParameters);
end

function generateMockData_normallyDistributed(bias, tID, nStudents, groupsA, groupsC, groupsN, tScores, fileString, pretestParameters, defaultImprovementParameters)

if bias > 1
    error('WARNING: Code not adapted to deal with biases larger than 1.')
end

out = [];

minPreTest = pretestParameters.min;
maxPreTest = pretestParameters.max;
avgPreTest = pretestParameters.avg;
stdPreTest = pretestParameters.std;

min_defaultImprovement = defaultImprovementParameters.min;
max_defaultImprovement = defaultImprovementParameters.max;
avg_defaultImprovement = defaultImprovementParameters.avg;
std_defaultImprovement = defaultImprovementParameters.std;

stored_defaultImprovement = [];
stored_correction = [];
stored_actualImprovement = [];
stored_groupsA = [];
stored_groupsC = [];
stored_groupsN = [];

for i = 1:length(tID)
    
    pre = genSemiRandom(minPreTest, avgPreTest, stdPreTest, nStudents(i), maxPreTest);
    defaultImprovement = genSemiRandom(min_defaultImprovement, avg_defaultImprovement, std_defaultImprovement, nStudents(i), max_defaultImprovement);
    
    if groupsA(i) == 1
        % A_top
        % addition between  0   and 2   for bias == 1, improvement of 1 on average
        % addition between -0.5 and 1.5 for bias == 0.5, improvement of 0.5 on average
        % addition between -0.9 and 1.1 for bias == 0.1, improvement of 0.1 on average
        minCorrection = bias - 1;
        maxCorrection = bias + 1;
        
    elseif groupsA(i) == 2
        % A_mid
        % addition between -1 and 1 for all biases, an (additional) improvement of 0 on average
        minCorrection = -1;
        maxCorrection = 1;
        
    elseif groupsA(i) == 3
        % A_bot
        % addition between -2   and 0   for bias == 1, decrease of 1 on average
        % addition between -1.5 and 0.5 for bias == 0.5, decreaseof 0.5 on average
        % addition between -1.1 and 0.9 for bias == 0.1, decrease of 0.1 on average
        minCorrection = -1 * (bias + 1);
        maxCorrection = -1 * (bias - 1);
        
    else
        error('Something went wrong in allocating groups')
    end
    
    avgCorrection = (minCorrection + maxCorrection) / 2;
    % A quarter of the range
    stdCorrection = 0.25 * (maxCorrection - minCorrection) / 2;
    correction = genSemiRandom(minCorrection, avgCorrection, stdCorrection, nStudents(i), maxCorrection);
    
    post = pre + defaultImprovement + correction;
    
    for j = 1:nStudents(i)
        tmp2(j,:) = tScores(i,:);
        tmp3(j,:) = tID(i);
    end
    
    percImpr = (post - pre) ./ pre .* 100;
    out = [out; post' pre' tmp3 tmp2 percImpr'];
    
    stored_defaultImprovement(end+1: end+length(defaultImprovement)) = defaultImprovement;
    stored_correction(end+1: end+length(defaultImprovement)) = correction;
    stored_actualImprovement(end+1: end+length(defaultImprovement)) = post - pre;
    stored_groupsA(end+1: end+length(defaultImprovement)) = groupsA(i);
    stored_groupsC(end+1: end+length(defaultImprovement)) = groupsC(i);
    stored_groupsN(end+1: end+length(defaultImprovement)) = groupsN(i);
    
    clear tmp tmp2 tmp3 post pre correction1
    
end

% Create BoxPlots of stored values
fileString_tmp = [fileString '_boxplot_defaultImprovement_bias(' num2str(bias) ')'];
ystring = 'Default Improvement (absolute points)';
storedVals_Boxplots(stored_defaultImprovement, fileString_tmp, stored_groupsA, ystring)

fileString_tmp = [fileString '_boxplot_improvementCorrection_bias(' num2str(bias) ')'];
ystring = 'Improvement Correction (absolute points)';
storedVals_Boxplots(stored_correction, fileString_tmp, stored_groupsA, ystring)

fileString_tmp = [fileString '_boxplot_simulatedImprovement_bias(' num2str(bias) ')'];
ystring = 'Simulated Improvement (absolute points)';
storedVals_Boxplots(stored_actualImprovement, fileString_tmp, stored_groupsA, ystring)
% Also store the improvement for the other groupings.
storedVals_Boxplots_CN(stored_actualImprovement, fileString_tmp, stored_groupsC, ystring, 'C')
storedVals_Boxplots_CN(stored_actualImprovement, fileString_tmp, stored_groupsN, ystring, 'N')

% NB: Could add a pre-post test plot, but as there is no effect of pre test
% score, this is not (yet) necessary.

output.data = out;
output.dataLabel = {'Post','Pre','tID','A1','C1','N1','A2','C2','N2','A','C','N', 'percImpr'}; % Where A = affective, C = continuance, N = normative
disp(output.dataLabel)

csvwrite([fileString '_mockData_normallyDistributed_' num2str(bias) '.csv'],out,1,0)

end

function number = genSemiRandom(minimumNumber, avg, std, n, maximumNumber)

whilLim = 0;
number = minimumNumber - 1;
secondStatement = 1;

while any(number < minimumNumber) || secondStatement
    
    number = normrnd(avg,std,[1 n]);
    
    whilLim = whilLim + 1;
    if whilLim > 10000
        error(['WARNING: With the current minimumNumber (' num2str(minimumNumber) ') (and possibly maximumNumber), even after 1000 iterations, no random set was found. Change the <minimumNumber>, or <avg>, or <std>...'])
    end
    if nargin == 5
        secondStatement = any(number > maximumNumber);
    else
        secondStatement = 0;
    end
end

end

function studentNumberBoxplots(nStudents, fileString, groups)
%% OVERALL
h = figure();
boxplot(nStudents)
ylabel('Number of Students per Teacher')
xticks('')
axis([0.75 1.25 0 ceil(max(nStudents)/10)*10])
h.PaperUnits = 'inches';
h.PaperPosition = [0 0 2.5 4];
print( h, '-r300' ,'-dpng' ,['randomGenInfo\' fileString '_boxplot_nStudents.png']) % here you can specify filename extensions
close(h)

%% PER TEACHER RANK
h = figure();
boxplot(nStudents, groups)
ylabel('Number of Students per Teacher')
set(gca,'XTick',[1 2 3],'XTickLabels',{'A_top','A_mid','A_bot'})
axis([0.5 3.5 0 ceil(max(nStudents)/10)*10])
h.PaperUnits = 'inches';
h.PaperPosition = [0 0 4.25 4];
print( h, '-r300' ,'-dpng' ,['randomGenInfo\' fileString '_boxplot_nStudents_per_A_rank.png']) % here you can specify filename extensions
close(h)
end

function storedVals_Boxplots(stored_vals, fileString, stored_groups, ystring)
roundYTicksTo = 3;
ytickmin = floor(min(stored_vals)/roundYTicksTo)*roundYTicksTo;
ytickmax = ceil(max(stored_vals)/roundYTicksTo)*roundYTicksTo;
%% OVERALL
h = figure();
boxplot(stored_vals)
ylabel(ystring)
xticks('')
axis([0.75 1.25 ytickmin ytickmax])
h.PaperUnits = 'inches';
h.PaperPosition = [0 0 2.5 4];
print( h, '-r300' ,'-dpng' ,['randomGenInfo\' fileString '.png']) % here you can specify filename extensions
close(h)

%% PER TEACHER RANK
h = figure();
boxplot(stored_vals, stored_groups)
ylabel(ystring)
set(gca,'XTick',[1 2 3],'XTickLabels',{'A_top','A_mid','A_bot'})
axis([0.5 3.5 ytickmin ytickmax])
h.PaperUnits = 'inches';
h.PaperPosition = [0 0 4.25 4];
print( h, '-r300' ,'-dpng' ,['randomGenInfo\' fileString '_per_A_rank.png']) % here you can specify filename extensions
close(h)
end

function storedVals_Boxplots_CN(stored_vals, fileString, stored_groups, ystring,ord)
roundYTicksTo = 3;
ytickmin = floor(min(stored_vals)/roundYTicksTo)*roundYTicksTo;
ytickmax = ceil(max(stored_vals)/roundYTicksTo)*roundYTicksTo;
% %% OVERALL
% h = figure();
% boxplot(stored_vals)
% ylabel(ystring)
% xticks('')
% axis([0.75 1.25 ytickmin ytickmax])
% h.PaperUnits = 'inches';
% h.PaperPosition = [0 0 2.5 4];
% print( h, '-r300' ,'-dpng' ,['randomGenInfo\' fileString '.png']) % here you can specify filename extensions
% close(h)

%% PER TEACHER RANK
h = figure();
boxplot(stored_vals, stored_groups)
ylabel(ystring)
set(gca,'XTick',[1 2 3],'XTickLabels',{[ord '_top'],[ord '_mid'],[ord '_bot']})
axis([0.5 3.5 ytickmin ytickmax])
h.PaperUnits = 'inches';
h.PaperPosition = [0 0 4.25 4];
print( h, '-r300' ,'-dpng' ,['randomGenInfo\' fileString '_per_' ord '_rank.png']) % here you can specify filename extensions
close(h)
end

function groups = getGroupScores(col, tScores)
[~,tmp] = sort(tScores(:,col));
steps = floor(length(tScores)/3);

groups = tmp;
groups(tmp(1:steps)) = 3; % bottom
groups(tmp(steps +1:steps *2)) = 2; % middle
groups(tmp(steps*2 +1:end)) = 1; % top

end

function tScores = createRandomTeacherScores(nTeachers)

tScores = NaN(nTeachers, 9);
for i = 1:6    
    tScores(:,i) = randi(35, nTeachers, 1);    
end
for i = 1:3
    tScores(:,i + 6) = sum(tScores(:,[i i+3]),2);
end
end
