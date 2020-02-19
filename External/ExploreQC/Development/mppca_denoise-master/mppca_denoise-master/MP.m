function [Xdn, sigma, npars] = MP(X, nbins, centering)

    % "MP": matrix denoising and noiseestimation by exploiting  data redundancy in the PCA domain using universal properties of the eigenspectrum of
    % random covariance matrices, i.e. Marchenko Pastur distribution
    %
    %  [Xdn, Sigma, npars] = MP(X, nbins)
    %       output:
    %           - Xdn: [MxN] denoised data matrix
    %           - sigma: [1x1] noise level
    %           - npars: [1x1] number of significant components
    %       input:
    %           - X: [MxN] data matrix
    %           - nbins: number of histogram bins for visualization. If
    %           empty or not provided, no graphs will be shown. 
    %
    %  Author: Jelle Veraart (jelle.veraart@nyumc.org)
    %  Copyright (c) 2016 New York University
    %       
    %      Permission is hereby granted, free of charge, to any non-commercial entity
    %      ('Recipient') obtaining a copy of this software and associated
    %      documentation files (the 'Software'), to the Software solely for
    %      non-commercial research, including the rights to use, copy and modify the
    %      Software, subject to the following conditions: 
    %       
    %        1. The above copyright notice and this permission notice shall be
    %      included by Recipient in all copies or substantial portions of the
    %      Software. 
    %       
    %        2. THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND,
    %      EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIESOF
    %      MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN
    %      NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BELIABLE FOR ANY CLAIM,
    %      DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
    %      OTHERWISE, ARISING FROM, OUT OF ORIN CONNECTION WITH THE SOFTWARE OR THE
    %      USE OR OTHER DEALINGS IN THE SOFTWARE. 
    %       
    %        3. In no event shall NYU be liable for direct, indirect, special,
    %      incidental or consequential damages in connection with the Software.
    %      Recipient will defend, indemnify and hold NYU harmless from any claims or
    %      liability resulting from the use of the Software by recipient. 
    %       
    %        4. Neither anything contained herein nor the delivery of the Software to
    %      recipient shall be deemed to grant the Recipient any right or licenses
    %      under any patents or patent application owned by NYU. 
    %       
    %        5. The Software may only be used for non-commercial research and may not
    %      be used for clinical care. 
    %       
    %        6. Any publication by Recipient of research involving the Software shall
    %      cite the references listed below.
    % 
    % REFERENCES
    %      Veraart, J.; Fieremans, E. & Novikov, D.S. Diffusion MRI noise mapping
    %      using random matrix theory Magn. Res. Med., 2016, early view, doi:
    %      10.1002/mrm.26059
    %      Veraart, J.; Novikov, D.S.; Christiaens, D.; Ades-Aron, B.; Sijbers, J. & Fieremans, E. 
    %      Denoising of diffusion MRI using random matrix theory, NeuroImage, Magn. Res. Med., 2016, early view, 
    %      DOI: 10.1016/j.neuroimage.2016.08.016
    

    
    if ~exist('nbins', 'var') || isempty(nbins)
        nbins=0;
    end
  
    [M, N] = size(X);
    
      
     
    if ~exist('centering', 'var') || isempty(centering)
        centering = false;
    end
    
    if centering
            colmean = mean(X, 1);
            X = X - repmat(colmean, [M, 1]);
    end 
            
    R = min(M, N);
    scaling = ones(R-centering, 1);
    if M>N
       %scaling = M/N;
       scaling = (M - (0:R-centering-1)) / N;
       scaling(scaling<1) = 1;
       scaling = scaling(:);
    end

    [u, vals, v] = svd(X, 'econ');
    vals = diag(vals).^2 / N;   
    csum = cumsum(vals(R-centering:-1:1)); cmean = csum(R-centering:-1:1)./(R-centering:-1:1)'; sigmasq_1 = cmean./scaling;
    
    gamma = (M - (0:R-centering-1)) / N;
    rangeMP = 4*sqrt(gamma(:));
    rangeData = vals(1:R-centering) - vals(R-centering);
    sigmasq_2 = rangeData./rangeMP;
    
    t = find(sigmasq_2 < sigmasq_1, 1);
    sigma = sqrt(sigmasq_1(t));
 
    npars = t-1; 
    
    if nbins>0

            [~, range] = MarchenkoPasturDistribution(rand(), sigma, M-npars, N);
            [p, ~] = MarchenkoPasturDistribution([range(1):diff(range)/100:range(2)], sigma, M-npars, N);
        
            figure;  
            hold on
            range_ = [vals(R-centering), vals(npars+1)];
            binwidth = diff(range_)/nbins;          % Finds the width of each bin.
            scale = M * binwidth;
            
            x = histc(vals(1:R-centering), [range_(1):diff(range_)/(nbins-1):range_(2)]);
            
            bar([range_(1):diff(range_)/(nbins-1):range_(2)], x/nansum(p))
            plot([range(1):diff(range)/100:range(2)], real(p)*scale/nansum(p), 'r', 'LineWidth', 3)
            
            
            xlabel('$\lambda$', 'FontName', 'Times', 'FontSize', 20, 'Interpreter', 'Latex')
            ylabel('$p(\lambda$)', 'FontName', 'Times', 'FontSize', 20, 'Interpreter', 'Latex')
            set(gca, 'FontSize', 20, 'box', 'on', 'LineWidth', 2, 'FontSize', 20);
            
            title(['sigma = ', num2str(sigma), ' and npars = ', num2str(npars)])
    end
    
    vals(t:R) = 0;
    Xdn = u*diag(sqrt(N*vals))*v';
    
    if centering
        Xdn = Xdn + repmat(colmean, [M, 1]);
    end
            

end

function [p, range] = MarchenkoPasturDistribution(lambda, sigma, M, N)
        Q = M/N;
        lambda_p = sigma^2*(1 + sqrt(Q)).^2;
        lambda_m = sigma^2*(1 - sqrt(Q)).^2;
        p = sqrt((lambda_p - lambda).*(lambda-lambda_m))./(2*pi*Q*lambda*sigma.^2);         
        p(lambda < lambda_m) = 0;       
        p(lambda > lambda_p) = 0;
        range = [lambda_m, lambda_p];
end