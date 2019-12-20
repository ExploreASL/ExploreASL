function [sexOut] = RandomizeSex(SampleSize,FemalePercentage)
% RandomizeSex Simple quick & dirty function to create pseudorandom
% sample with M/F distribution per samplesize & percentage of Females
% 
    FemaleCount = round(FemalePercentage*SampleSize); % count number of females
    sexT(1:FemaleCount,1) = 1; % create part of row of females
    sexT(FemaleCount+1:SampleSize,1) = 0; % rest of the row is males
    sexT(:,2) = rand([1 SampleSize])'; % pseudocreate random order
    sexT = sortrows(sexT,2); % sort according to pseudorandom order
    sexOut = sexT(:,1)'; % output column
end