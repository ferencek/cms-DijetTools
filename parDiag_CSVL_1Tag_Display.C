{
  TMatrixDSym M = TMatrixDSym(4);

  double cov_matrix[] = {
      1.8533e-13,  2.5379e-08, -2.3285e-09,  9.2719e-10,
      2.5379e-08,   0.0056911,  0.00010015,  0.00024742,
     -2.3285e-09,  0.00010015,  0.00020257,  5.9752e-05,
      9.2719e-10,  0.00024742,  5.9752e-05,  3.6641e-05
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
  double fit_pars[] = { 1.54555e-05, 7.86212e+00, 5.40672e+00, 4.97943e-02 };
  p.SetMatrixArray(fit_pars);
  cout << "Fit parameters in original basis:" << endl;
  p.Print();

  TMatrixD p_diag = T_tr*p;
  cout << "Fit parameters in diagonal basis:" << endl;
  p_diag.Print();
}