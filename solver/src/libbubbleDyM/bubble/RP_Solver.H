/* solving the Rayleigh-Plesset */
/* ======================================================= */

    
    const label n = dt/1e-8;
    const scalar endTime = dt;
     // Create the ODE system
    RPODE ode(pc,p0,rhoc,muc,pv,sigma,R0_);
    dictionary dict;
    dict.add("solver", "RKF45");
    // Create the selected ODE system solver
    autoPtr<ODESolver> odeSolver = ODESolver::New(ode, dict);
    // Initialise the ODE system fields
    scalar xStart = mesh().time().value();
 
    scalar dx;
    // time step
    if (n == 0)
    {
        dx = 0;
    }
    else
    {
        dx = endTime/n;
    }
    // Initialise displacement and velocity
    scalarField yStart(ode.nEqns());
    yStart[0]= R_;
    yStart[1]= dRt_;
    
    //integration initial step
    scalar dxEst = 1e-16;
    scalar xEnd = 0;

    // Print the evolution of the solution and the time-step
    scalarField dyStart(ode.nEqns());

    //Integration loop
    for(label i=0; i<n;i++)
    {
         xEnd = xStart+dx;
         ode.derivatives(xStart,yStart,dyStart);
	  
         odeSolver->solve(xStart,xEnd,yStart,dxEst);
	     xStart = xEnd;
         if(mag(yStart[0] - 1e-8) < 1e-8)
         {
               yStart[0] = 1e-8;
               yStart[1] = 0;
               break; 
         }


    }
    R_ = yStart[0];
    dRt_ = yStart[1];
    d_ = 2.0*R_;
    scalar kappa = 4.0/3.0;

     if (R_ >= R0_)
     {
          kappa=1;
     }
    scalar pB_ = (p0 +2*sigma/R0_-pv)*Foam::pow(R0_/R_,3*kappa)+pv;
    ddRt_ = -3*pow(dRt_,2)/(2*R_) + (pB_- pc)/(rhoc*R_) - (2*sigma + 4*muc*dRt_)/(rhoc*sqr(R_));
    if(mag(R_ - 1e-8) < 1e-8)
    {
        td.keepParticle = false;
        R_ = 1e-8;
        dRt_ = 0.0;
        d_ = 2.0*R_;
        ddRt_ = 0.0;
    }

/* ======================================================= */    
