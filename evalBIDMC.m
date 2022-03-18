clear all,clc,close all

load('bidmc_data.mat')

absERR = zeros(size(data,2),8);
for i = 1:size(data,2)
    for j = 1:8
        ECG = data(i).ekg.v((j-1)*60*125+1:j*60*125);
        PPG = data(i).ppg.v((j-1)*60*125+1:j*60*125);
        fs = data(i).ekg.fs;
        breathREF = mean(data(i).ref.params.rr.v((j-1)*60+1:j*60),'omitnan');
        breathVAL = RRestimate(fs,PPG,ECG);
        absERR(i,j) = abs(breathREF-breathVAL);
    end
end
MAE = mean(absERR,'all','omitnan');
figure
bar(mean(absERR,2,'omitnan'))
hold on
line([1 53],[MAE MAE])
title(['MAE: ' num2str(MAE) '; median: ' num2str(median(reshape(absERR,[],1)))])
figure
boxplot(absERR')
figure
boxplot(reshape(absERR,[],1))
