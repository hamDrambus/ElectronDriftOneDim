Double_t fit_f_1(Double_t *x,Double_t *par) {
    if (x[0]>par[0])
      return par[1]*TMath::Exp(-(x[0]-par[0])/par[2]);
    return 0;
}

Double_t fit_f_2(Double_t *x,Double_t *par) {
    if (x[0]>par[0])
      return par[1]*TMath::Exp(-(x[0]-par[0])/par[2]) + par[3]*TMath::Exp(-(x[0]-par[0])/par[4]);
    return 0;
}
