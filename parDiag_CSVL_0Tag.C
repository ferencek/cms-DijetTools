{
  TMatrixDSym M = TMatrixDSym(4);

  double cov_matrix[] = {
      1.468e-05,  1.365e-04, -1.254e-05,  4.920e-06,
      1.365e-04,  2.095e-03,  3.625e-05,  8.898e-05,
     -1.254e-05,  3.625e-05,  7.328e-05,  2.148e-05,
      4.920e-06,  8.898e-05,  2.148e-05,  1.315e-05
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
  double fit_pars[] = { 2.28255e-01, 8.51130e+00, 5.42168e+00,  5.23767e-02 };
  p.SetMatrixArray(fit_pars);
  cout << "Fit parameters in original basis:" << endl;
  p.Print();

  TMatrixD p_diag = T_tr*p;
  cout << "Fit parameters in diagonal basis:" << endl;
  p_diag.Print();
}