%% Vectorized MATLAB implementation of the STAPLE algorithm by Warfield et al.
% This code currently only supports the case of binary segmentations.
% 
% Function:  [W, p, q] = STAPLE(D)
% Parameter: D, data matrix of segmentations, dimensions VOXELS x RATERS
% Returns:   W, est. weight matrix for each voxel
%            p, est. vector of sensitivities for each expert
%            q, est. vector of specificities for each expert
%
% You can simply threshold the resulting matrix W to get an estimated ground truth segmentation,
% e.g. T = (W >= .5)
%
% Literature: Warfield, Simon K., Kelly H. Zou, and William M. Wells. 
%            "Simultaneous truth and performance level estimation (STAPLE): 
%            an algorithm for the validation of image segmentation." 
%           Medical Imaging, IEEE Transactions on 23.7 (2004): 903-921.
%
% Andreas Husch
% 2013-07-10, 2015-07-06, 2016-05-10
% mail@andreashusch.de
% Mohammad Haft-Javaherian (mh973@cornell.edu) 
% Fall 2017
%% Example for usage:
% %Using test data (2D prostate segmentations) from original publication
% s1 = nrrdread('seg001.nrrd');
% s2 = nrrdread('seg002.nrrd');
% s3 = nrrdread('seg003.nrrd');
% s4 = nrrdread('seg004.nrrd');
% s5 = nrrdread('seg005.nrrd');
% imageDims = size(s1);
% D = [s1(:), s2(:), s3(:), s4(:), s5(:)]; % pixels in rows, raters in columns
% [W, p, q]= STAPLE(D);
% % p,q values of your raters:
% p
% q
% Estimated ground truth image:
% gtImage = reshape((W >= .5), imageDims);
% figure, imagesc(gtImage)


%% Vectorized implementation of the classical STAPLE-Algorithm by Warfield et al. for binary segmentations
function [W] = STAPLE_FWD(D, p, q)
    %% Inputs & Checks

     originalSize=size(D);
     vessel_ID = find(any(~isnan(D),2));
     user_ID = find(any(~isnan(D),1));
     D=D(vessel_ID, user_ID);
     D = double(D);
     p = p(user_ID)';
     q = q(user_ID)';
     N = size(D,2); %Number of raters
     
    %% Parameters

    
    % Initial sensitivity and specificity parameter p(j),q(j) for all
    % raters j
    DnotNaN = ~ismissing(D);

    Tprior = (sum(D(:),'omitnan')/sum(DnotNaN(:))); % note dependence on (sub)volume size, final result depends on this prior (which is not an implementation issue but a basic limitation of the EM approach)

    W = ones(1,size(D,1));
    [row,col] = size(D);
    
    
    %% EM
 
        % The following  code is equivalent to this loop by MUCH faster
        %     for i = 1:length(D)
        %         W(i) = ((prod(p(D(i,:))) * prod(1 - p(~D(i,:)))) * Tprior) / ... 
        %                ((prod(p(D(i,:))) * prod(1 - p(~D(i,:)))) * Tprior + (prod(q(~D(i,:))) * prod(1 - q(D(i,:))))) * (1- Tprior) ;
        %         %NOTE that prod([]) = 1
        %     end       
        P = repmat(p,size(D,1), 1);
        Q = repmat(q,size(D,1), 1);
        P_given_D = P .* D; %TODO: use bsxfun instead of repmat?
        P_given_D(P_given_D(:) == 0) = 1; %
        Q_given_D = 1 - Q .* D;
        Q_given_D(Q_given_D(:) == 0) = 1; % alternative: initialise with 1 and set Q_given_D(D) = 1- P(D) 
        tempDP = D;
        tempDP(DnotNaN) = P(DnotNaN).* ~D(DnotNaN);
        compP_given_not_D  = 1 - tempDP;
        compP_given_not_D(compP_given_not_D(:)== 0) = 1;
        tempDQ = D;
        tempDQ(DnotNaN) = Q(DnotNaN).* ~D(DnotNaN);
        compQ_given_not_D  = tempDQ;
        compQ_given_not_D(compQ_given_not_D(:)== 0) = 1;
        %hist(prod(P_given_D','omitnan'),100)
        % W(i) can be interpreted as the prob. of voxel i being true (i.e. is part of the ground-truth y) for given p(1:N), q(1:N) 
        prodP_given_D = prod(P_given_D','omitnan');
        prodP_given_D(prodP_given_D == 0) = 10^-200;
        prodCompP = prod(compP_given_not_D','omitnan');
        prodCompP(prodCompP == 0) = 10^-200;
        prodQ_given_D = prod(Q_given_D','omitnan');
        prodQ_given_D(prodQ_given_D == 0) = 10^-200;
        prodCompQ = prod(compQ_given_not_D','omitnan');
        prodCompQ(prodCompQ == 0) = 10^-200;
        a = prodP_given_D .* prodCompP;
        a(a == 0) = 10^-200;
        b = prodQ_given_D .* prodCompQ;
        b(b == 0) = 10^-200;
        W = (a * Tprior) ./ ...
           ((a * Tprior) + b * (1 - Tprior)); %#ok<UDIM>
        W0 = W;
        W = nan(originalSize(1), 1);
        W(vessel_ID) = W0;
end