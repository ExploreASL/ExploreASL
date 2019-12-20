    ACA_L                  = uint8( vasc(:,:,:,1)==1;
    ACA_R                  = uint8( vasc(:,:,:,2)==1;
    
    MCA_L                  = uint8( vasc(:,:,:,3)==1;
    MCA_R                  = uint8( vasc(:,:,:,4)==1;
    
    PCA_L                  = uint8( vasc(:,:,:,5)==1;
    PCA_R                  = uint8( vasc(:,:,:,6)==1;
    
    %% Randomize
    
    RandomN                 = rand(500,1);
    
    for iR=1:80
        % select random image
        IndexN(iR,1)              = round( size(ASL_untreated.Data.data,1) * RandomN(iR)/(max(RandomN)-min(RandomN)) );
    end
    
    IndexN                  = unique(IndexN);
    
    
    for iR=1:40
        
        
        