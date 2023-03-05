%% Output Params

set(0,'DefaultFigureColor',APSslideColor)

t = datetime('now');
t.Format = "yyMMdd";
dat_str0 = string(t);

if saveIntResults
    % save all trials in a single directory at highest level
    tryName = num2str(numTrials);
    trial_dir = strcat(dat_str0,'_',tryName,'trialRuns');
    mkdir(trial_dir)
else
    % no group directory required
    trial_dir = '..';
end
%%

%results folder for this particular data run

t = datetime('now');
t.Format = "yyMMddhhmm";
dat_str = string(t);
dir_name = strcat(trial_dir,'/',dat_str,'_',num2str(L(1)),'spins3D');
mkdir(dir_name)
mkdir(dir_name,'frames')
mkdir(dir_name,'png')
mkdir(dir_name,'fig')
mkdir(dir_name,'txt')
