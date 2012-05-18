{
  TMatrixDSym M = TMatrixDSym(4);

  double cov_matrix[] = {
      6.9288e-16,  4.9494e-09, -4.7562e-10,  1.7933e-10,
      4.9494e-09,    0.057985,  0.00093532,   0.0025374,
     -4.7562e-10,  0.00093532,   0.0021434,  0.00062798,
      1.7933e-10,   0.0025374,  0.00062798,  0.00038357
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
  double fit_pars[] = { 2.91837e-07, 5.46300e+00, 6.37524e+00, 2.15882e-01 };
  p.SetMatrixArray(fit_pars);
  cout << "Fit parameters in original basis:" << endl;
  p.Print();

  TMatrixD p_diag = T_tr*p;
  cout << "Fit parameters in diagonal basis:" << endl;
  p_diag.Print();
}