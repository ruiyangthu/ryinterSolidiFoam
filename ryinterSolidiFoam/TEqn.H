{
    int iter = 0;
    scalar residual = 1;
//    scalar meanResidual = 1;

   do
    {
       iter++;
 
        alphal.storePrevIter();
 
        cp = alpha1*(alphal*cpl + alphas*cps) + alpha2*cpg;
        lambda = alpha1*(alphal*lambdal + alphas*lambdas) + alpha2*lambdag;
        phiCp = fvc::interpolate(cp)*phi;


        fvScalarMatrix pcTEqn
        (
                fvm::ddt(cp, T)
          + fvm::div(phiCp, T)
          - fvm::laplacian(lambda/rho, T)
          
         + alpha1*L*fvc::ddt(alphal)
          + alpha1*L*fvc::div(phi, alphal)
          - alpha1*Tmelt*fvc::ddt(cp)
          - alpha1*Tmelt*fvc::div(phiCp)
 
        );

        pcTEqn.solve();


        Tcorr = (Tl-Ts)*alphal+Ts;  
        alphal = max(min(alphal+alphaRel*cp/L*(T-Tcorr),scalar(1)),scalar(0));
        alphas = max(min(scalar(1) - alphal ,scalar(1)),scalar(0));
        residual = max(mag(alphal.internalField()-alphal.prevIter().internalField()));
       // meanResidual = sum(mag(alphal.internalField()-alphal.prevIter().internalField())*mesh.V())/sum(mesh.V()).value();
          
        Info<< "Liquid Melt volume fraction alphal = "
        << alphal.weightedAverage(mesh.Vsc()).value()
        << "  Min(alphal) = " << min(alphal).value()
        << "  Max(alphal) = " << max(alphal).value()
        << endl;

        Info<< "Solid Melt volume fraction alphas = "
        << alphas.weightedAverage(mesh.Vsc()).value()
        << "  Min(alphas) = " << min(alphas).value()
        << "  Max(alphas) = " << max(alphas).value()
        << "  iter    = " << iter
        << endl;
    } while ((iter < minTCorr || residual > alphaTol ) && iter <= maxTCorr); 
//while ((residual > alphaTol ) );

    alphas = max(min(alphas,scalar(1)),scalar(0));
    h = cp*(T-Tmelt)+alphal*alpha1*L;
    DC = DCl*pow(alphas,2)/(pow(alphal,3)+DCs);
    
//    rhok = scalar(1)-max(beta*(T-Tl),scalar(0));
}
