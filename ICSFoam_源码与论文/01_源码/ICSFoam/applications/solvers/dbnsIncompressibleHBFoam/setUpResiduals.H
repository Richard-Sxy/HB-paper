// One scalar equation (continuity) and one vector equation (momentum).
// The global pressure/velocity blocks both use phase index K.
residualsIO initResI(1, 1);

initResI.sInitRes()[0] = initRes().sInitRes()[K];
initResI.vInitRes()[0] = initRes().vInitRes()[K];

initResI.sFinalRes()[0] = initRes().sFinalRes()[K];
initResI.vFinalRes()[0] = initRes().vFinalRes()[K];

solnControl[K].setResidual(initResI);
