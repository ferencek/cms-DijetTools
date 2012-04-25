{
  TMatrixDSym M = TMatrixDSym(4);

  double cov_matrix[] = {
      1.537e-05,  2.288e-04, -2.130e-05,  8.337e-06,
      2.288e-04,  5.569e-03,  9.516e-05,  2.434e-04,
     -2.130e-05,  9.516e-05,  2.011e-04,  5.917e-05,
      8.337e-06,  2.434e-04,  5.917e-05,  3.619e-05
  };

  M.SetMatrixArray(cov_matrix);
  cout << "Covariance matrix:" << endl;
  M.Print();
  cout << "Determinant: " << M.Determinant() << endl << endl;
  
  const TMatrixDSymEigen eigen(M);
  TVectorD e_val = eigen.GetEigenValues();
  cout << "Covariance matrix eigenvalues:" << endl;
  e_val.Print();
  cout << "Errors in diagonal basis:" << endl;
  e_val.Sqrt();
  e_val.Print();
  TMatrixD T = eigen.GetEigenVectors();
  cout << "Transformation matrix T (formed from the covariance matrix eigenvectors):" << endl;
  T.Print();

  TMatrixD T_tr = T;
  T_tr.T();
  cout << "Transpose of T:" << endl;
  T_tr.Print();

  TMatrixD p = TMatrixD(4,1);
  double fit_pars[] = { 1.41453e-01, 8.47691e+00, 4.97184e+00, -3.60922e-02 };
  p.SetMatrixArray(fit_pars);
  cout << "Fit parameters in original basis:" << endl;
  p.Print();

  TMatrixD p_diag = T_tr*p;
  cout << "Fit parameters in diagonal basis:" << endl;
  p_diag.Print();
}