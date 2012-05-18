{
  TMatrixDSym M = TMatrixDSym(4);

  double cov_matrix[] = {
      6.2459e-13,  2.8414e-08, -2.5629e-09,  1.0336e-09,
      2.8414e-08,   0.0021325,  3.8518e-05,  9.0851e-05,
     -2.5629e-09,  3.8518e-05,  7.3834e-05,  2.1758e-05,
      1.0336e-09,  9.0851e-05,  2.1758e-05,  1.3361e-05
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
  double fit_pars[] = { 4.68740e-05, 8.52546e+00, 5.40954e+00, 4.96618e-02 };
  p.SetMatrixArray(fit_pars);
  cout << "Fit parameters in original basis:" << endl;
  p.Print();

  TMatrixD p_diag = T_tr*p;
  cout << "Fit parameters in diagonal basis:" << endl;
  p_diag.Print();
}