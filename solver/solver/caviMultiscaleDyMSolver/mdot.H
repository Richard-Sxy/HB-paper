dimensionedScalar deltaT = mesh.time().deltaT();

//volScalarField mdotbeta = Cvb * rho1*rho2/rho * (Sa + fvc::div(phi,alpha1)-alpha1*fvc::div(phi)); 
volScalarField mdotbeta = (-Cvb) * rho1*rho2/rho * Sa ; 

volScalarField vdotbeta = mag(mdotbeta * (1.0/rho1 - 1.0/rho2))/(mag(p-pSat)+ 0.001*pSat);

//volScalarField vDotbetaP = vdotbeta*neg(beta1 - 1.0) + vDotvP * neg(alpha1 - 0.999);
volScalarField vDotbetaP =  vDotvP  * neg(alpha1 - 0.999);


