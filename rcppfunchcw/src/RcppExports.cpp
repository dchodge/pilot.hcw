// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// ab_kin_normal
void ab_kin_normal(NumericVector& predicted_titres, const NumericVector& theta, const List& infection_info, const List& vaccination_info, const List& setup_data, const List& indexing, const List& antigenic_maps, const List& other_pars);
RcppExport SEXP _rcppfunchcw_ab_kin_normal(SEXP predicted_titresSEXP, SEXP thetaSEXP, SEXP infection_infoSEXP, SEXP vaccination_infoSEXP, SEXP setup_dataSEXP, SEXP indexingSEXP, SEXP antigenic_mapsSEXP, SEXP other_parsSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector& >::type predicted_titres(predicted_titresSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const List& >::type infection_info(infection_infoSEXP);
    Rcpp::traits::input_parameter< const List& >::type vaccination_info(vaccination_infoSEXP);
    Rcpp::traits::input_parameter< const List& >::type setup_data(setup_dataSEXP);
    Rcpp::traits::input_parameter< const List& >::type indexing(indexingSEXP);
    Rcpp::traits::input_parameter< const List& >::type antigenic_maps(antigenic_mapsSEXP);
    Rcpp::traits::input_parameter< const List& >::type other_pars(other_parsSEXP);
    ab_kin_normal(predicted_titres, theta, infection_info, vaccination_info, setup_data, indexing, antigenic_maps, other_pars);
    return R_NilValue;
END_RCPP
}
// ab_kin_vac
void ab_kin_vac(NumericVector& predicted_titres, const NumericVector& theta, const List& infection_info, const List& vaccination_info, const List& setup_data, const List& indexing, const List& antigenic_maps, const List& other_pars);
RcppExport SEXP _rcppfunchcw_ab_kin_vac(SEXP predicted_titresSEXP, SEXP thetaSEXP, SEXP infection_infoSEXP, SEXP vaccination_infoSEXP, SEXP setup_dataSEXP, SEXP indexingSEXP, SEXP antigenic_mapsSEXP, SEXP other_parsSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector& >::type predicted_titres(predicted_titresSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const List& >::type infection_info(infection_infoSEXP);
    Rcpp::traits::input_parameter< const List& >::type vaccination_info(vaccination_infoSEXP);
    Rcpp::traits::input_parameter< const List& >::type setup_data(setup_dataSEXP);
    Rcpp::traits::input_parameter< const List& >::type indexing(indexingSEXP);
    Rcpp::traits::input_parameter< const List& >::type antigenic_maps(antigenic_mapsSEXP);
    Rcpp::traits::input_parameter< const List& >::type other_pars(other_parsSEXP);
    ab_kin_vac(predicted_titres, theta, infection_info, vaccination_info, setup_data, indexing, antigenic_maps, other_pars);
    return R_NilValue;
END_RCPP
}
// ab_kin_vac_general
void ab_kin_vac_general(NumericVector& predicted_titres, const NumericVector& theta, const List& infection_info, const List& vaccination_info, const List& setup_data, const List& indexing, const List& antigenic_maps, const List& other_pars);
RcppExport SEXP _rcppfunchcw_ab_kin_vac_general(SEXP predicted_titresSEXP, SEXP thetaSEXP, SEXP infection_infoSEXP, SEXP vaccination_infoSEXP, SEXP setup_dataSEXP, SEXP indexingSEXP, SEXP antigenic_mapsSEXP, SEXP other_parsSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector& >::type predicted_titres(predicted_titresSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const List& >::type infection_info(infection_infoSEXP);
    Rcpp::traits::input_parameter< const List& >::type vaccination_info(vaccination_infoSEXP);
    Rcpp::traits::input_parameter< const List& >::type setup_data(setup_dataSEXP);
    Rcpp::traits::input_parameter< const List& >::type indexing(indexingSEXP);
    Rcpp::traits::input_parameter< const List& >::type antigenic_maps(antigenic_mapsSEXP);
    Rcpp::traits::input_parameter< const List& >::type other_pars(other_parsSEXP);
    ab_kin_vac_general(predicted_titres, theta, infection_info, vaccination_info, setup_data, indexing, antigenic_maps, other_pars);
    return R_NilValue;
END_RCPP
}
// ab_kin_vac_prev_hist
void ab_kin_vac_prev_hist(NumericVector& predicted_titres, const NumericVector& theta, const List& infection_info, const List& vaccination_info, const List& setup_data, const List& indexing, const List& antigenic_maps, const List& other_pars);
RcppExport SEXP _rcppfunchcw_ab_kin_vac_prev_hist(SEXP predicted_titresSEXP, SEXP thetaSEXP, SEXP infection_infoSEXP, SEXP vaccination_infoSEXP, SEXP setup_dataSEXP, SEXP indexingSEXP, SEXP antigenic_mapsSEXP, SEXP other_parsSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector& >::type predicted_titres(predicted_titresSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const List& >::type infection_info(infection_infoSEXP);
    Rcpp::traits::input_parameter< const List& >::type vaccination_info(vaccination_infoSEXP);
    Rcpp::traits::input_parameter< const List& >::type setup_data(setup_dataSEXP);
    Rcpp::traits::input_parameter< const List& >::type indexing(indexingSEXP);
    Rcpp::traits::input_parameter< const List& >::type antigenic_maps(antigenic_mapsSEXP);
    Rcpp::traits::input_parameter< const List& >::type other_pars(other_parsSEXP);
    ab_kin_vac_prev_hist(predicted_titres, theta, infection_info, vaccination_info, setup_data, indexing, antigenic_maps, other_pars);
    return R_NilValue;
END_RCPP
}
// ab_kin_vac_log_normal
void ab_kin_vac_log_normal(NumericVector& predicted_titres, const NumericVector& theta, const List& infection_info, const List& vaccination_info, const List& setup_data, const List& indexing, const List& antigenic_maps, const List& other_pars);
RcppExport SEXP _rcppfunchcw_ab_kin_vac_log_normal(SEXP predicted_titresSEXP, SEXP thetaSEXP, SEXP infection_infoSEXP, SEXP vaccination_infoSEXP, SEXP setup_dataSEXP, SEXP indexingSEXP, SEXP antigenic_mapsSEXP, SEXP other_parsSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector& >::type predicted_titres(predicted_titresSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const List& >::type infection_info(infection_infoSEXP);
    Rcpp::traits::input_parameter< const List& >::type vaccination_info(vaccination_infoSEXP);
    Rcpp::traits::input_parameter< const List& >::type setup_data(setup_dataSEXP);
    Rcpp::traits::input_parameter< const List& >::type indexing(indexingSEXP);
    Rcpp::traits::input_parameter< const List& >::type antigenic_maps(antigenic_mapsSEXP);
    Rcpp::traits::input_parameter< const List& >::type other_pars(other_parsSEXP);
    ab_kin_vac_log_normal(predicted_titres, theta, infection_info, vaccination_info, setup_data, indexing, antigenic_maps, other_pars);
    return R_NilValue;
END_RCPP
}
// ab_kin_vac_set_point
void ab_kin_vac_set_point(NumericVector& predicted_titres, const NumericVector& theta, const List& infection_info, const List& vaccination_info, const List& setup_data, const List& indexing, const List& antigenic_maps, const List& other_pars);
RcppExport SEXP _rcppfunchcw_ab_kin_vac_set_point(SEXP predicted_titresSEXP, SEXP thetaSEXP, SEXP infection_infoSEXP, SEXP vaccination_infoSEXP, SEXP setup_dataSEXP, SEXP indexingSEXP, SEXP antigenic_mapsSEXP, SEXP other_parsSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector& >::type predicted_titres(predicted_titresSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const List& >::type infection_info(infection_infoSEXP);
    Rcpp::traits::input_parameter< const List& >::type vaccination_info(vaccination_infoSEXP);
    Rcpp::traits::input_parameter< const List& >::type setup_data(setup_dataSEXP);
    Rcpp::traits::input_parameter< const List& >::type indexing(indexingSEXP);
    Rcpp::traits::input_parameter< const List& >::type antigenic_maps(antigenic_mapsSEXP);
    Rcpp::traits::input_parameter< const List& >::type other_pars(other_parsSEXP);
    ab_kin_vac_set_point(predicted_titres, theta, infection_info, vaccination_info, setup_data, indexing, antigenic_maps, other_pars);
    return R_NilValue;
END_RCPP
}
// ab_kin_normal_rel
void ab_kin_normal_rel(NumericVector& predicted_titres, const NumericVector& theta, const List& infection_info, const List& vaccination_info, const List& setup_data, const List& indexing, const List& antigenic_maps, const List& other_pars);
RcppExport SEXP _rcppfunchcw_ab_kin_normal_rel(SEXP predicted_titresSEXP, SEXP thetaSEXP, SEXP infection_infoSEXP, SEXP vaccination_infoSEXP, SEXP setup_dataSEXP, SEXP indexingSEXP, SEXP antigenic_mapsSEXP, SEXP other_parsSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector& >::type predicted_titres(predicted_titresSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const List& >::type infection_info(infection_infoSEXP);
    Rcpp::traits::input_parameter< const List& >::type vaccination_info(vaccination_infoSEXP);
    Rcpp::traits::input_parameter< const List& >::type setup_data(setup_dataSEXP);
    Rcpp::traits::input_parameter< const List& >::type indexing(indexingSEXP);
    Rcpp::traits::input_parameter< const List& >::type antigenic_maps(antigenic_mapsSEXP);
    Rcpp::traits::input_parameter< const List& >::type other_pars(other_parsSEXP);
    ab_kin_normal_rel(predicted_titres, theta, infection_info, vaccination_info, setup_data, indexing, antigenic_maps, other_pars);
    return R_NilValue;
END_RCPP
}
// ab_kin_vac_rel
void ab_kin_vac_rel(NumericVector& predicted_titres, const NumericVector& theta, const List& infection_info, const List& vaccination_info, const List& setup_data, const List& indexing, const List& antigenic_maps, const List& other_pars);
RcppExport SEXP _rcppfunchcw_ab_kin_vac_rel(SEXP predicted_titresSEXP, SEXP thetaSEXP, SEXP infection_infoSEXP, SEXP vaccination_infoSEXP, SEXP setup_dataSEXP, SEXP indexingSEXP, SEXP antigenic_mapsSEXP, SEXP other_parsSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector& >::type predicted_titres(predicted_titresSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const List& >::type infection_info(infection_infoSEXP);
    Rcpp::traits::input_parameter< const List& >::type vaccination_info(vaccination_infoSEXP);
    Rcpp::traits::input_parameter< const List& >::type setup_data(setup_dataSEXP);
    Rcpp::traits::input_parameter< const List& >::type indexing(indexingSEXP);
    Rcpp::traits::input_parameter< const List& >::type antigenic_maps(antigenic_mapsSEXP);
    Rcpp::traits::input_parameter< const List& >::type other_pars(other_parsSEXP);
    ab_kin_vac_rel(predicted_titres, theta, infection_info, vaccination_info, setup_data, indexing, antigenic_maps, other_pars);
    return R_NilValue;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_rcppfunchcw_ab_kin_normal", (DL_FUNC) &_rcppfunchcw_ab_kin_normal, 8},
    {"_rcppfunchcw_ab_kin_vac", (DL_FUNC) &_rcppfunchcw_ab_kin_vac, 8},
    {"_rcppfunchcw_ab_kin_vac_general", (DL_FUNC) &_rcppfunchcw_ab_kin_vac_general, 8},
    {"_rcppfunchcw_ab_kin_vac_prev_hist", (DL_FUNC) &_rcppfunchcw_ab_kin_vac_prev_hist, 8},
    {"_rcppfunchcw_ab_kin_vac_log_normal", (DL_FUNC) &_rcppfunchcw_ab_kin_vac_log_normal, 8},
    {"_rcppfunchcw_ab_kin_vac_set_point", (DL_FUNC) &_rcppfunchcw_ab_kin_vac_set_point, 8},
    {"_rcppfunchcw_ab_kin_normal_rel", (DL_FUNC) &_rcppfunchcw_ab_kin_normal_rel, 8},
    {"_rcppfunchcw_ab_kin_vac_rel", (DL_FUNC) &_rcppfunchcw_ab_kin_vac_rel, 8},
    {NULL, NULL, 0}
};

RcppExport void R_init_rcppfunchcw(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
