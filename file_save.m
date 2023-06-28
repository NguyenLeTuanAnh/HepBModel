% save results as separate files
% get the directory of results files
pathname = fileparts('/home/pca/Documents/Group A Strep/GAS_transmission_model/GAS_Skin_Throat_model/results/');
[stamp] = timestamp;
filebase = strcat('/GAS_',num2str(scenario.num_iter),'i_s', ...
    num2str(scenario.setting), stamp);
filesummary = strcat(pathname, filebase, '_summary.csv');
filetheta = strcat(pathname, filebase,'_theta.csv');
filell = strcat(pathname, filebase,'_ll.csv');
csvwrite(filesummary, results.summaries)
csvwrite(filetheta, results.theta)
csvwrite(filell, results.ll)


