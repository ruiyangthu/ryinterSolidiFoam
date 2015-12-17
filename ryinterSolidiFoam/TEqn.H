{
    int iter = 0;
     scalar residual = 1;
    // scalar meanResidual = 1;

//     volScalarField fl ("fl", max(min((T-Ts)/(Tl-Ts),scalar(1)),scalar(0)));
  //   volScalarField fs ("fs", scalar(1) - fl);    
//predict
     cp = alpha1*fl*cpl + alpha1*fs*cps + alpha2*cpg;
     lambda = alpha1*fl*lambdal + alpha1*fs*lambdas + alpha2*lambdag;
     solve (fvm::ddt(cp, T)+ fvm::div(phiCp, T)- fvm::laplacian(lambda/rho, T));
     Info<< "Teq predict "  << endl;
     fl = max(min((T-Ts)/(Tl-Ts),scalar(1)),scalar(0));
     fs = max(min(scalar(1) - fl ,scalar(1)),scalar(0));
   
//correction
    do
    {
        iter++;
        fl.storePrevIter();

       //         cp = fl*cpl + fs*cps + alpha2*cpg;
        cp = alpha1*fl*cpl + alpha1*fs*cps + alpha2*cpg;
        phiCp = fvc::interpolate(cp)*phi;
        lambda = alpha1*fl*lambdal + alpha1*fs*lambdas + alpha2*lambdag;

        phiCp = fvc::interpolate(cp)*phi;

        fvScalarMatrix TEqn
        (
            fvm::ddt(cp, T)
          + fvm::div(phiCp, T)
          - fvm::laplacian(lambda/rho, T)
          ==
            L*fvc::ddt(alpha1,fl)
          + L*fvc::div(phi, alpha1*fl)
          - Tmelt*fvc::ddt(cp)
          - Tmelt*fvc::div(phiCp)
          
        );

        TEqn.solve();

        Tcorr = (Tl-Ts)*fl+Ts;  
        fl = max(min(fl+alphaRel*cp/L*(T-Tcorr),scalar(1)),scalar(0));
        fs = max(min(scalar(1) - fl ,scalar(1)),scalar(0));
        residual = max(mag(fl.internalField()-fl.prevIter().internalField()));
       // meanResidual = sum(mag(fl.internalField()-fl.prevIter().internalField())*mesh.V())/sum(mesh.V()).value();
          
        Info<< "Liquid Melt volume fraction fl = "
        << fl.weightedAverage(mesh.Vsc()).value()
        << "  Min(fl) = " << min(fl).value()
        << "  Max(fl) = " << max(fl).value()
        << endl;

        Info<< "Solid Melt volume fraction fs = "
        << fs.weightedAverage(mesh.Vsc()).value()
        << "  Min(fs) = " << min(fs).value()
        << "  Max(fs) = " << max(fs).value()
        << "  iter    = " << iter
        << endl;
    } while ((iter < minTCorr || residual > alphaTol ) && iter <= maxTCorr);
    alphas = max(min(fs,scalar(1)),scalar(0));
    h = cp*(T-Tmelt)+fl*alpha1*L;
    DC = DCl*pow(fs,2)/(pow(fl,3)+DCs);
    
//    rhok = scalar(1)-max(beta*(T-Tl),scalar(0));
}
