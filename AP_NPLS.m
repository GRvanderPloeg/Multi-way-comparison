% PARAFAC functionality
addpath("../N-way-shell/Scripts/"); % own scripts
addpath("../N-way-shell/N-way toolbox/"); % from Rasmus Bro

%%
% Load raw cytokines data
cytokines_data = readmatrix("./Data/AP/input_deduplicated_RvdP.csv", Filetype="delimitedtext", Delimiter=" ");
cytokines_meta = readmatrix("./Data/AP/input_deduplicated_metadata_RvdP.csv", Filetype="delimitedtext", Delimiter=" ", OutputType="string");

onlyCases = cytokines_meta(:,6) == "case";
cytokines_data = cytokines_data(onlyCases,:);
cytokines_meta = cytokines_meta(onlyCases,:);

path = "./Models/NPLS/AP/";

%%
% Attach pain/no-pain metadata to these case subjects
meta2 = readmatrix("./Data/AP/Root_meta_data_parafac.txt", Filetype="delimitedtext", Delimiter="\t", OutputType="string");
meta2 = meta2(:,[2,8]);

cytokines_meta(:,7) = "s";
allIndividuals = unique(cytokines_meta(:,1));

for i=1:length(allIndividuals)
    subjectID = allIndividuals(i);
    status = meta2(meta2(:,1) == subjectID, 2);
    if length(status) > 1 % somehow there are duplicate entries
        status = status(1);
    end
    cytokines_meta(cytokines_meta(:,1) == subjectID, 7) = status;
end

%%
% Log transform
[I,J] = size(cytokines_data);
vectorizedX = reshape(cytokines_data, I*J, 1);
pseudocount = min(vectorizedX(vectorizedX>0));

cytokines_data_log = log(cytokines_data + pseudocount);

%%
% Reshape to 3-way matrix
keepIndividuals = true;
cytokines_cube = rawDataToCube(cytokines_data_log, cytokines_meta(:,1), cytokines_meta(:,2), keepIndividuals);

%%
% Center and scale
[cytokines_cnt, cytokines_means] = centerData(cytokines_cube, 1);
[cytokines_cnt_scl, cytokines_stds] = scaleData(cytokines_cnt, 2);

%%
% Might want to add editing the metadata here
[I,J,K] = size(cytokines_cnt_scl);

subjectMeta = cytokines_meta(:,[1 7]);
subjectMeta = sortrows(subjectMeta, 1);
subjectMeta = unique(subjectMeta, "rows", "stable"); % stable stops resorting rows

subjectMeta(subjectMeta(:,2) == "A", 2) = "Asymptomatic";
subjectMeta(subjectMeta(:,2) == "S", 2) = "Symptomatic";

%featureMeta = [1:J; ones(1,J,1)]';
featureMeta = ["VEGF" "CRP" "GM-CSF" "IL1alpha" "IL1beta" "IL4" "IL6" "IL8" "IL10" "IL12p70" "IL17A" "IFNgamma" "MIP1alpha" "OPG" "TNFalpha" "RANKL"]';

%%
% Initialize options
Options = [];
Options(1) = 1e-6;  % convergence (default 1e-6)
Options(2) = 2;     % initialization (default 1)
Options(3) = 0;     % resulting plot (default 0)
const = [0 0 0];

%%
% Filter out outliers
df_cnt_scl_filtered = cytokines_cnt_scl;
subjectMeta_filtered = subjectMeta;
featureMeta_filtered = featureMeta;

timeMeta =  ["Before extraction" "Before extraction" "Before extraction" "After extraction" "After extraction", "After extraction"]';
timeMeta_filtered = timeMeta;

% Remove outlier subjects
% df_cnt_scl_filtered([25],:,:) = [];
% subjectMeta_filtered([25],:) = [];

% Remove timepoint 6
%df_cnt_scl_filtered(:,:,6) = [];

%% prepare y

Y = zeros(size(subjectMeta,1), 1);
Y(subjectMeta(:,2) == "Asymptomatic",:) = 1;
Y = Y - mean(Y);

%%
% Determine correct number of NPLS components
maxComponents = 10;
[XValResult_saliva,~] = ncrossreg(cytokines_cnt_scl, Y, maxComponents, 0);

%%
% Plot RMSEP of CV and save it
%set(gcf, 'Units', 'Normalized', 'outerposition', [0, 0, 1, 1], 'PaperUnits', 'Inches', 'PaperSize', [7.25, 9.125])
%plot(1:maxComponents, XValResult.RMSEP); xlim([1 maxComponents]); title(""); xlabel("Num. components"); ylabel("RMSEP"); yyaxis right; plot(1:maxComponents, XValResult.Percent.Xexp(:,1)); ylim([0 100]); ylabel("Variance explained (%)");
%saveas(gcf, path+"RMSEP_CV.jpg");
%close();
XValResult = XValResult_saliva;

% NEW APPROACH
df = [XValResult.RMSEP' XValResult.Percent.Xexp(:,1) XValResult.Percent.Yexp(:,1) XValResult.PRESS'];
% 
% plot(1:maxComponents, XValResult.RMSEP); xlim([1 maxComponents]); title(""); xlabel("Num. components"); ylabel("RMSEP"); yyaxis right; plot(1:maxComponents, XValResult.Percent.Xexp(:,1)); ylim([0 100]); ylabel("Variance explained (%)");
% plot(1:maxComponents, XValResult.Percent.Xexp(:,1)); xlim([1 maxComponents]); title(""); xlabel("Num. components"); ylabel("Variance explained (%)"); plot(1:maxComponents, XValResult.Percent.Yexp(:,1));
% 
% colororder({'r', 'b', 'm'});
%plot(1:maxComponents, df(:,[2 3])); xlim([1 maxComponents]); title(""); xlabel("Num. components"); ylabel("Variance explained (%)"); yyaxis right; plot(1:maxComponents, df(:,1)); ylabel("RMSEP"); legend("VarExpX", "VarExpY", "RMSEP");

% ALTERNATIVE
plot(1:maxComponents, df(:,[2 3])); xlim([1 maxComponents]); title(""); xlabel("Num. components"); ylabel("Variance explained (%)"); legend("X", "Y");
saveas(gcf, path+"RMSEP_CV_new_varExps.jpg");
close();
plot(1:maxComponents, df(:,1)); xlim([1 maxComponents]); title(""); xlabel("Num. components"); ylabel("RMSEP");
saveas(gcf, path+"RMSEP_CV_new_RMSEP.jpg");
close();
plot(1:maxComponents, df(:,4)); xlim([1 maxComponents]); title(""); xlabel("Num. components"); ylabel("PRESS");
saveas(gcf, path+"RMSEP_CV_new_PRESS.jpg");
close();

%%
% Run NPLS model
numComponents = 1;
[Xfactors,Yfactors,Core,B,ypred,ssx,ssy,reg] = npls(cytokines_cnt_scl, Y, numComponents);

%%
% Save NPLS models
metaData = {subjectMeta_filtered, featureMeta_filtered, timeMeta_filtered};
annotatedModel = annotateModel(cytokines_cnt_scl, Xfactors, metaData);
savePARAFAC(cytokines_cnt_scl, Xfactors, annotatedModel, path + "AP");

%%
% Save ypred of NPLS models
%y_saliva = [ypred_saliva saliva_subjectMeta_filtered testosterone];
y = [ypred(:,:,1) subjectMeta_filtered Y];
writematrix(y, path + "AP_" + numComponents + "_ypred.csv");

%%
% Save coefficients of NPLS models
for i=1:size(reg, 2)
    writematrix(reg{i}, path + "AP_" + numComponents + "_coeff_" + i + ".csv");
end

%%
% Save crossvalidated coefficients of the NPLS models
[coeff_mean, coeff_std] = CV_coeff_NPLS(cytokines_cnt_scl, Y, numComponents, 1);
writematrix(coeff_mean, path + "AP_" + numComponents + "_CVcoeff_means.csv");
writematrix(coeff_std, path + "AP_" + numComponents + "_CVcoeff_stds.csv");