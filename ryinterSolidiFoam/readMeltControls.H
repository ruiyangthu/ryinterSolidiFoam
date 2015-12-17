const dictionary& meltingDict = mesh.solutionDict().subDict("MELTING");
int minTCorr(readLabel(meltingDict.lookup("minTempCorrector")));
int maxTCorr(readLabel(meltingDict.lookup("maxTempCorrector")));
scalar alphaTol(readScalar(meltingDict.lookup("alphaTolerance")));
scalar alphaRel(readScalar(meltingDict.lookup("alphaRelaxatio")));
