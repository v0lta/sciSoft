% Make the plots from the article.
% About
%   Author       - Peter Opsomer (Peter.Opsomer@cs.kuleuven.be)
%   History      - Created October 2013, last edit February 2015
%% Initialising
format longe; close all; clear variables;
set(0,'DefaultFigureWindowStyle','docked');

for fig = 1:2
    switch fig
        case 1
            alpha = 0; beta = 0; h = @(x) exp(-7*x.^(2*2) );
			maxOrder = 7;maxP2 = 8;
        case 2
            alpha = -1/2; beta = 0; h = @(x) 1./sqrt(x+3); 
			maxOrder = 6;maxP2 = 7;
    end
    tic;
    [P,gammaP,alphaP,betaP] = exactPolys(alpha,beta,h,2^maxP2);
    timeOPQ = toc
    
    %% Computing higher order terms
    tic;
    [c, d, Dinf, psi, dpsi] = contour_integrals(alpha,beta,h,maxOrder);
    [Uright,Uleft] = UQ(alpha,beta,Dinf,c,d,maxOrder);
    timePrecompute = toc
    
    %% Computing exact results and asymptotic expansions
    zInt = 0.2+0.5*1i;
    zL = -0.97;
    zOut = zInt; % Should lie within contour !
    
    piEvi = zeros(maxP2,1); piEvbL = zeros(maxP2,1); piEvo = zeros(maxP2,1);
	ipiEvi = zeros(maxP2,maxOrder); bpiLEv = zeros(maxP2,maxOrder); opiEv = zeros(maxP2,maxOrder);
	gammas = zeros(maxP2,1);
	for tn = 1:maxP2
		gammas(tn) = gammaP(2^tn+1)*2^(2^tn);
	end
	
    for tn = 1:maxP2
        piEvi(tn) = P(zInt,2^tn +1)/gammas(tn);
        piEvbL(tn) = P(zL,2^tn +1)/gammas(tn);
        piEvo(tn) = P(zOut,2^tn +1)/gammas(tn);
        for i = 1:maxOrder
            ipiEvi(tn,i) = asy_lens(2^tn,zInt,alpha,beta,h,psi,i,Dinf,Uright,Uleft);
            bpiLEv(tn,i) = asy_left(2^tn,zL,alpha,beta,h,psi,i,Dinf,Uright,Uleft);
            opiEv(tn,i) = asy_outer(2^tn,zOut,alpha,beta,h,psi,i,Dinf,Uright,Uleft);
        end
    end
    
    %% Testing order of convergence
    switch fig
        case 1
            plotConv(piEvbL,bpiLEv,'Left disk');
        case 2
            plotConv(piEvi,ipiEvi,'Lens');
            plotConv(piEvo,opiEv,'Outside region');
    end
    shg
end
set(0, 'DefaultFigureVisible', 'on');
