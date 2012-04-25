{
  TMatrixDSym M = TMatrixDSym(4);

  double cov_matrix[] = {
      1.329e-07,  6.671e-05, -6.534e-06,  2.428e-06,
      6.671e-05,  5.443e-02,  8.487e-04,  2.441e-03,
     -6.534e-06,  8.487e-04,  2.082e-03,  6.104e-04,
      2.428e-06,  2.441e-03,  6.104e-04,  3.715e-04
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
  double fit_pars[] = { 4.11375e-03, 6.37762e+00, 5.58514e+00,  4.93530e-02 };
  p.SetMatrixArray(fit_pars);
  cout << "Fit parameters in original basis:" << endl;
  p.Print();

  TMatrixD p_diag = T_tr*p;
  cout << "Fit parameters in diagonal basis:" << endl;
  p_diag.Print();
}