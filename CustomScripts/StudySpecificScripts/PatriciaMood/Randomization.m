for iS=1:13
    for iSess=1:19
        Randomization{(iS-1)*19+iSess,1}     = Sex{iS,1};
        Randomization{(iS-1)*19+iSess,2}     = ['ASL_' num2str(iSess)];
    end
end

