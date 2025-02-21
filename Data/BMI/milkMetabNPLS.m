% PARAFAC functionality
addpath("./N-way-shell\Scripts\"); % own scripts
addpath("./N-way-shell\N-way toolbox/"); % from Rasmus Bro

%%
% Load data
df = readmatrix("./milkMetabNumeric.csv");
taxonomy = readmatrix("./milk_metab_CAS_numbers.csv", OutputType="string");
taxonomy = taxonomy(2:end,:);
sampleInfo = readmatrix("./milkMetab_sampleMeta.csv", Filetype="text", OutputType="string");

sampleInfo2 = readmatrix("./milkMetab_sampleInfo_w_Anthropometrics.csv", Filetype="text", OutputType="string");
%sampleInfo2 = readmatrix("./20240906_data/milk_metab_meta.csv", OutputType="string");
% all(sampleInfo(:,1) == sampleInfo2(:,1)) samples are ordered the same

% Filter taxa to taxa with at least 1 non-zero value
featureMask = (sum(df) > 0);
df = df(:, featureMask);
taxonomy = taxonomy(featureMask,:);

% Filter subjects to be present in X and Y
% The column order is as follows for the data from 20240906
% 1 is RCID
% 2 is ppBMI
% 3 is BMI.group
% 4 is BMI.delta.3m
% 12 is delta.whz
% 13 is BMIz.1y
% 14 is BMIz.2y
% 15 is BMIz.3y

% For the ordinary data it is
% 2 is BMI
% 14 is WHZ.6m
selection = 14;

mask = sampleInfo2(:,selection) ~= "NA";
df = df(mask,:);
sampleInfo = sampleInfo(mask,:);
sampleInfo2 = sampleInfo2(mask,:);
Y = unique([sampleInfo(:,11) sampleInfo2(:,selection)], 'rows');
Y = strrep(Y, ",", "."); % replace comma with point as decimal separator

path = "./20241107_milkMetab_NPLS_BMI_2comp/";
%%
% Log transform
df_log = log(df);

%%
% Make into cube
keepIndividuals = true;
[df_cube, subjectMeta, conditionMeta] = rawDataToCube(df_log, sampleInfo(:,11), str2double(sampleInfo(:,4)), keepIndividuals);

%%
% Center and scale
df_cnt = centerData(df_cube, 1);
df_cnt_scl = scaleData(df_cnt, 2);

%%
% Manually remove subject 56 
% Prepare metadata and remove outliers
subjectMeta_filtered = unique(sampleInfo(:,[11 3]), "rows");
featureMeta_filtered = taxonomy;
timeMeta_filtered = conditionMeta;

df_cnt_scl_filtered = df_cnt_scl;

% Prepare Y
Ynum = str2double(Y(:,2)); % confirmed to be ordered the same way as subjectMeta_filtered
Ycnt = Ynum - mean(Ynum);

%%
% Initialize options
Options = [];
Options(1) = 1e-6;  % convergence (default 1e-6)
Options(2) = 2;     % initialization (default 1)
Options(3) = 0;     % resulting plot (default 0)
const = [0 0 0];

%%
% Determine correct number of NPLS components per niche
maxComponents = 10;
[XValResult,~] = ncrossreg(df_cnt_scl_filtered, Ycnt, maxComponents, 0);

%%
% Plot RMSEP of CV and save it
%set(gcf, 'Units', 'Normalized', 'outerposition', [0, 0, 1, 1], 'PaperUnits', 'Inches', 'PaperSize', [7.25, 9.125])
%plot(1:maxComponents, XValResult.RMSEP); xlim([1 maxComponents]); title(""); xlabel("Num. components"); ylabel("RMSEP"); yyaxis right; plot(1:maxComponents, XValResult.Percent.Xexp(:,1)); ylim([0 100]); ylabel("Variance explained (%)");
%saveas(gcf, path+"RMSEP_CV.jpg");
%close();

% NEW APPROACH
df = [XValResult.RMSEP' XValResult.Percent.Xexp(:,1) XValResult.Percent.Yexp(:,1) XValResult.PRESS'];

%plot(1:maxComponents, XValResult.RMSEP); xlim([1 maxComponents]); title(""); xlabel("Num. components"); ylabel("RMSEP"); yyaxis right; plot(1:maxComponents, XValResult.Percent.Xexp(:,1)); ylim([0 100]); ylabel("Variance explained (%)");
%plot(1:maxComponents, XValResult.Percent.Xexp(:,1)); xlim([1 maxComponents]); title(""); xlabel("Num. components"); ylabel("Variance explained (%)"); plot(1:maxComponents, XValResult.Percent.Yexp(:,1));

%colororder({'r', 'b', 'm'});
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
% Find optimal number of components based on RMSEP CV
%[bestRMSEP_saliva, numComponents_saliva] = min(XValResult_saliva.RMSEP);

% OVERRIDE
numComponents = 1;
bestRMSEP = XValResult.RMSEP(numComponents);

%%
% Run NPLS
[Xfactors,Yfactors,Core,B,ypred,ssx,ssy,reg] = npls(df_cnt_scl_filtered, Ycnt, numComponents);

%%
% Save NPLS models
metaData = {subjectMeta_filtered, featureMeta_filtered, timeMeta_filtered};
annotatedModel = annotateModel(df_cnt_scl_filtered, Xfactors, metaData);
savePARAFAC(df_cnt_scl_filtered, Xfactors, annotatedModel, path + "Rasmus");

%%
% Save ypred of NPLS models
%y_saliva = [ypred_saliva saliva_subjectMeta_filtered testosterone];
y = [ypred(:,:,1) subjectMeta_filtered Ycnt];
writematrix(y, path + "Rasmus_" + numComponents + "_ypred.csv");

%%
% Save coefficients of NPLS models
for i=1:size(reg, 2)
    writematrix(reg{i}, path + "Rasmus_" + numComponents + "_coeff_" + i + ".csv");
end

%%
% Save crossvalidated coefficients of the NPLS models
[coeff_mean, coeff_std] = CV_coeff_NPLS(df_cnt_scl_filtered, Ycnt, numComponents, 1);
writematrix(coeff_mean, path + "Rasmus_" + numComponents + "_CVcoeff_means.csv");
writematrix(coeff_mean, path + "Rasmus_" + numComponents + "_CVcoeff_stds.csv");