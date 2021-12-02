
#include <Rcpp.h>
using namespace Rcpp;

#include <Rcpp.h>
using namespace Rcpp;

#ifndef MAX
#define MAX(a,b) ((a) < (b) ? (b) : (a)) // define MAX function for use later
#endif

#ifndef MIN
#define MIN(a,b) ((a) < (b) ? (a) : (b)) // define MAX function for use later
#endif


//' Base boosting fast
//' 
//' A fast implementation of the basic boosting function, giving predicted titres for a number of samples for one individual. Note that this version attempts to minimise memory allocations.
//' @family boosting_functions
//' @seealso \code{\link{titre_data_fast}}
// [[Rcpp::export]]
void ab_kin_normal(
            NumericVector &predicted_titres,
            const NumericVector &theta,
            const List &infection_info,
            const List &vaccination_info,
            const List &setup_data,
            const List &indexing,
            const List &antigenic_maps,
            const List &other_pars
          ){

  int max_vaccinations;
  bool vac_flag = !vaccination_info["vac_null_ind"];
  NumericVector vaccination_times;
  if (vac_flag) {
    NumericVector vaccination_times = vaccination_info["vac_times"];
    IntegerVector vaccination_strain_indices_tmp = vaccination_info["vac_indices"];
    max_vaccinations = vaccination_times.size();
  }

  NumericVector infection_times = infection_info["inf_times"];
  IntegerVector infection_strain_indices_tmp = infection_info["inf_indices"];

  int index_in_samples = indexing["index_in_samples"];
  int end_index_in_samples = indexing["end_index_in_samples"];
  int start_index_in_data = indexing["start_index_in_data"];

  NumericVector sample_times = setup_data["sample_times"];
  IntegerVector measurement_strain_indices = setup_data["measurement_strain_indices"];
  IntegerVector nrows_per_blood_sample = setup_data["nrows_per_blood_sample"];
  int number_strains = setup_data["number_strains"];

  double sampling_time;
  double time;
  double n_inf;
  double n_vac;
  int x_inf;
  int x_vac;

  double wane_amount, wane_amount_vac;
  double seniority;
  double vac_suppress;

  int n_titres;

  int end_index_in_data;
  int tmp_titre_index;
  int inf_map_index;
  int vac_map_index;
  int index;

  int max_infections = infection_times.size();

  double mu = theta["mu"];
  double mu_short = theta["mu_short"];
  double wane = theta["wane"];
  double tau = theta["tau"];
  double kappa = theta["kappa"];

  NumericVector antigenic_map_long = antigenic_maps["long"];
  NumericVector antigenic_map_short = antigenic_maps["short"];

  NumericVector inf_vac_times;
  if (vac_flag) {
    if (max_vaccinations > 0  & max_infections > 0) {
      inf_vac_times = union_(infection_times, vaccination_times);
    } else if(max_vaccinations == 0  & max_infections > 0) {
      inf_vac_times = infection_times;
    } else if (max_vaccinations > 0  & max_infections == 0) {
      inf_vac_times = vaccination_times;
    } 
  }
  else {
    inf_vac_times = infection_times;
  }
  std::sort(inf_vac_times.begin(), inf_vac_times.end());

  int max_inf_vac_times = inf_vac_times.size();

  LogicalVector indicatior_inf(max_inf_vac_times);
  LogicalVector indicatior_vac(max_inf_vac_times);

  for (int i = 0; i < inf_vac_times.size(); i++) {
    if (vac_flag) {
      for (int k = 0; k < vaccination_times.size(); k++) {
        indicatior_vac[i] = inf_vac_times[i] == vaccination_times[k];
        if (indicatior_vac[i])
          break;
      }
    }
    for (int j = 0; j < infection_times.size(); j++) {
      indicatior_inf[i] = inf_vac_times[i] == infection_times[j];
      if (indicatior_inf[i])
        break;
    }
  }

  for(int j = index_in_samples; j <= end_index_in_samples; ++j){
    sampling_time = sample_times[j];
    n_inf = 0.0;
    n_vac = 0.0;
    x_inf = 0.0;
    x_vac = 0.0;
    // Find number of titres in the predicted_titres vector that correspond to this sample
    n_titres = nrows_per_blood_sample[j];
    // Only iterate through indices for this sample
    end_index_in_data = start_index_in_data + n_titres;
    tmp_titre_index = start_index_in_data;

    // Sum all infections that would contribute towards observed titres at this time
    for (int x = 0; x < max_inf_vac_times; ++x){
      if (indicatior_inf[x]) {
        ++n_inf;
        // sampling through each predicted infection
        // Only go further if this sample happened after the infection
        if(sampling_time > infection_times[x_inf]) {
            time = sampling_time - infection_times[x_inf]; // Time between sample and infection
            wane_amount = MAX(0, 1.0 - (wane*time)); // Basic waning function
            seniority = MAX(0, 1.0 - tau*(n_inf + n_vac - 1.0)); // Antigenic seniority
            
            inf_map_index = infection_strain_indices_tmp[x_inf]; // Index of this infecting strain in antigenic map
            for(int k = 0; k < n_titres; ++k){
              index = measurement_strain_indices[tmp_titre_index + k]*number_strains + inf_map_index;
              predicted_titres[tmp_titre_index + k] += seniority *
                ((mu*antigenic_map_long[index]) + (mu_short*antigenic_map_short[index])*wane_amount);
            }
        }
        ++x_inf;
      }
    }
    start_index_in_data = end_index_in_data;
  }
}

//' Base boosting fast
//' 
//' A fast implementation of the basic boosting function, giving predicted titres for a number of samples for one individual. Note that this version attempts to minimise memory allocations.
//' @family boosting_functions
//' @seealso \code{\link{titre_data_fast}}
// [[Rcpp::export]]
void ab_kin_vac(
            NumericVector &predicted_titres,
            const NumericVector &theta,
            const List &infection_info,
            const List &vaccination_info,
            const List &setup_data,
            const List &indexing,
            const List &antigenic_maps,
            const List &other_pars
          ){

  int max_vaccinations;
  bool vac_flag = !vaccination_info["vac_null_ind"];
  NumericVector vaccination_times;
  IntegerVector vaccination_strain_indices_tmp;
  if (vac_flag) {
    vaccination_times = vaccination_info["vac_times"];
    vaccination_strain_indices_tmp = vaccination_info["vac_indices"];
    max_vaccinations = vaccination_times.size();
  }

  NumericVector infection_times = infection_info["inf_times"];
  IntegerVector infection_strain_indices_tmp = infection_info["inf_indices"];

  int index_in_samples = indexing["index_in_samples"];
  int end_index_in_samples = indexing["end_index_in_samples"];
  int start_index_in_data = indexing["start_index_in_data"];

  NumericVector sample_times = setup_data["sample_times"];
  IntegerVector measurement_strain_indices = setup_data["measurement_strain_indices"];
  IntegerVector nrows_per_blood_sample = setup_data["nrows_per_blood_sample"];
  int number_strains = setup_data["number_strains"];

  double sampling_time;
  double time;
  double n_inf;
  double n_vac;
  int x_inf;
  int x_vac;

  double wane_amount, wane_amount_vac;
  double seniority;
  double vac_suppress;

  int n_titres;

  int end_index_in_data;
  int tmp_titre_index;
  int inf_map_index;
  int vac_map_index;
  int index;

  int max_infections = infection_times.size();

  double mu = theta["mu"];
  double mu_short = theta["mu_short"];
  double wane = theta["wane"];
  double tau = theta["tau"];
  double kappa = theta["kappa"];
  double mu_vac = theta["mu_vac"];
  double tau_prev_vac = theta["tau_prev_vac"];
  double mu_short_vac = theta["mu_short_vac"];
  double wane_vac = theta["wane_vac"];

  NumericVector antigenic_map_long = antigenic_maps["long"];
  NumericVector antigenic_map_short = antigenic_maps["short"];
  NumericVector antigenic_map_long_vac = antigenic_maps["long_vac"];
  NumericVector antigenic_map_short_vac = antigenic_maps["short_vac"];

  NumericVector inf_vac_times;
  if (vac_flag) {
    if (max_vaccinations > 0  & max_infections > 0) {
      inf_vac_times = union_(infection_times, vaccination_times);
    } else if(max_vaccinations == 0  & max_infections > 0) {
      inf_vac_times = infection_times;
    } else if (max_vaccinations > 0  & max_infections == 0) {
      inf_vac_times = vaccination_times;
    } 
  }
  else {
    inf_vac_times = infection_times;
  }
  std::sort(inf_vac_times.begin(), inf_vac_times.end());

  int max_inf_vac_times = inf_vac_times.size();

  LogicalVector indicatior_inf(max_inf_vac_times);
  LogicalVector indicatior_vac(max_inf_vac_times);

  for (int i = 0; i < inf_vac_times.size(); i++) {
    if (vac_flag) {
      for (int k = 0; k < vaccination_times.size(); k++) {
        indicatior_vac[i] = inf_vac_times[i] == vaccination_times[k];
        if (indicatior_vac[i])
          break;
      }
    }
    for (int j = 0; j < infection_times.size(); j++) {
      indicatior_inf[i] = inf_vac_times[i] == infection_times[j];
      if (indicatior_inf[i])
        break;
    }
  }

  for(int j = index_in_samples; j <= end_index_in_samples; ++j){
    sampling_time = sample_times[j];
    n_inf = 0.0;
    n_vac = 0.0;
    x_inf = 0.0;
    x_vac = 0.0;
    // Find number of titres in the predicted_titres vector that correspond to this sample
    n_titres = nrows_per_blood_sample[j];
    // Only iterate through indices for this sample
    end_index_in_data = start_index_in_data + n_titres;
    tmp_titre_index = start_index_in_data;

    // Sum all infections that would contribute towards observed titres at this time
    for (int x = 0; x < max_inf_vac_times; ++x){
      if (indicatior_inf[x]) {
        ++n_inf;
        // sampling through each predicted infection
        // Only go further if this sample happened after the infection
        if(sampling_time > infection_times[x_inf]) {
            time = sampling_time - infection_times[x_inf]; // Time between sample and infection
            wane_amount = MAX(0, 1.0 - (wane*time)); // Basic waning function
            seniority = MAX(0, 1.0 - tau*(n_inf + n_vac - 1.0)); // Antigenic seniority
            
            inf_map_index = infection_strain_indices_tmp[x_inf]; // Index of this infecting strain in antigenic map
            for(int k = 0; k < n_titres; ++k){
              index = measurement_strain_indices[tmp_titre_index + k]*number_strains + inf_map_index;
              predicted_titres[tmp_titre_index + k] += seniority *
                ((mu*antigenic_map_long[index]) + (mu_short*antigenic_map_short[index])*wane_amount);
            }
        }
        ++x_inf;
      }

      if (indicatior_vac[x] & (vac_flag)) {
          ++n_vac;
          if(sampling_time > vaccination_times[x_vac]) {
            time = sampling_time - vaccination_times[x_vac]; // Time er vaccination
            wane_amount_vac = MAX(0, 1.0 - (wane_vac*time)); // Basic waning function
            seniority = MAX(0, 1.0 - tau_prev_vac*(n_inf + n_vac - 1.0)); // Antigenic seniority
            vac_map_index = vaccination_strain_indices_tmp[x_vac]; // Index of this vaccinating strain in antigenic map

            for(int k = 0; k < n_titres; ++k){
              index = measurement_strain_indices[tmp_titre_index + k]*number_strains + vac_map_index;
              predicted_titres[tmp_titre_index + k] += (seniority) * ((mu*mu_vac * antigenic_map_long_vac[index]) +
                (mu_short_vac * antigenic_map_short_vac[index]) * wane_amount_vac);
            }
          }
          ++x_vac;
        }
    }
    start_index_in_data = end_index_in_data;
  }
}


//' Base boosting fast
//' 
//' A fast implementation of the basic boosting function, giving predicted titres for a number of samples for one individual. Note that this version attempts to minimise memory allocations.
//' @family boosting_functions
//' @seealso \code{\link{titre_data_fast}}
// [[Rcpp::export]]
void ab_kin_vac_prev_hist(
            NumericVector &predicted_titres,
            const NumericVector &theta,
            const List &infection_info,
            const List &vaccination_info,
            const List &setup_data,
            const List &indexing,
            const List &antigenic_maps,
            const List &other_pars
          ){

  int max_vaccinations;
  bool vac_flag = !vaccination_info["vac_null_ind"];
  NumericVector vaccination_times;
  IntegerVector vaccination_strain_indices_tmp;
  if (vac_flag) {
    vaccination_times = vaccination_info["vac_times"];
    vaccination_strain_indices_tmp = vaccination_info["vac_indices"];
    max_vaccinations = vaccination_times.size();
  }

  NumericVector infection_times = infection_info["inf_times"];
  IntegerVector infection_strain_indices_tmp = infection_info["inf_indices"];

  int index_in_samples = indexing["index_in_samples"];
  int end_index_in_samples = indexing["end_index_in_samples"];
  int start_index_in_data = indexing["start_index_in_data"];

  NumericVector sample_times = setup_data["sample_times"];
  IntegerVector measurement_strain_indices = setup_data["measurement_strain_indices"];
  IntegerVector nrows_per_blood_sample = setup_data["nrows_per_blood_sample"];
  int number_strains = setup_data["number_strains"];

  double sampling_time;
  double time;
  double n_inf;
  double n_vac;
  int x_inf;
  int x_vac;

  double wane_amount, wane_amount_vac;
  double seniority;
  double vac_suppress;

  int n_titres;

  int end_index_in_data;
  int tmp_titre_index;
  int inf_map_index;
  int vac_map_index;
  int index;

  int max_infections = infection_times.size();

  double mu = theta["mu"];
  double mu_short = theta["mu_short"];
  double wane = theta["wane"];
  double tau = theta["tau"];
  double kappa = theta["kappa"];
  double mu_vac = theta["mu_vac"];
  double tau_prev_vac = theta["tau_prev_vac"];
  double mu_short_vac = theta["mu_short_vac"];
  double wane_vac = theta["wane_vac"];
  double rho_boost = theta["rho_boost"];
  double rho_wane = theta["rho_wane"];

  double rho_boost_par;
  double rho_wane_par;

  NumericVector antigenic_map_long = antigenic_maps["long"];
  NumericVector antigenic_map_short = antigenic_maps["short"];
  NumericVector antigenic_map_long_vac = antigenic_maps["long_vac"];
  NumericVector antigenic_map_short_vac = antigenic_maps["short_vac"];

  NumericVector inf_vac_times;
  if (vac_flag) {
    if (max_vaccinations > 0  & max_infections > 0) {
      inf_vac_times = union_(infection_times, vaccination_times);
    } else if(max_vaccinations == 0  & max_infections > 0) {
      inf_vac_times = infection_times;
    } else if (max_vaccinations > 0  & max_infections == 0) {
      inf_vac_times = vaccination_times;
    } 
  }
  else {
    inf_vac_times = infection_times;
  }
  std::sort(inf_vac_times.begin(), inf_vac_times.end());

  int max_inf_vac_times = inf_vac_times.size();

  LogicalVector indicatior_inf(max_inf_vac_times);
  LogicalVector indicatior_vac(max_inf_vac_times);

  for (int i = 0; i < inf_vac_times.size(); i++) {
    if (vac_flag) {
      for (int k = 0; k < vaccination_times.size(); k++) {
        indicatior_vac[i] = inf_vac_times[i] == vaccination_times[k];
        if (indicatior_vac[i])
          break;
      }
    }
    for (int j = 0; j < infection_times.size(); j++) {
      indicatior_inf[i] = inf_vac_times[i] == infection_times[j];
      if (indicatior_inf[i])
        break;
    }
  }

  for(int j = index_in_samples; j <= end_index_in_samples; ++j){
    sampling_time = sample_times[j];
    n_inf = 0.0;
    n_vac = 0.0;
    x_inf = 0.0;
    x_vac = 0.0;
    // Find number of titres in the predicted_titres vector that correspond to this sample
    n_titres = nrows_per_blood_sample[j];
    // Only iterate through indices for this sample
    end_index_in_data = start_index_in_data + n_titres;
    tmp_titre_index = start_index_in_data;

    // Sum all infections that would contribute towards observed titres at this time
    for (int x = 0; x < max_inf_vac_times; ++x){
      if (indicatior_inf[x]) {
        ++n_inf;
        // sampling through each predicted infection
        // Only go further if this sample happened after the infection
        if(sampling_time > infection_times[x_inf]) {
            time = sampling_time - infection_times[x_inf]; // Time between sample and infection
            wane_amount = MAX(0, 1.0 - (wane*time)); // Basic waning function
            seniority = MAX(0, 1.0 - tau*(n_inf + n_vac - 1.0)); // Antigenic seniority
            
            inf_map_index = infection_strain_indices_tmp[x_inf]; // Index of this infecting strain in antigenic map
            for(int k = 0; k < n_titres; ++k){
              index = measurement_strain_indices[tmp_titre_index + k]*number_strains + inf_map_index;
              predicted_titres[tmp_titre_index + k] += seniority *
                ((mu*antigenic_map_long[index]) + (mu_short*antigenic_map_short[index])*wane_amount);
            }
        }
        ++x_inf;
      }

      if (indicatior_vac[x] & (vac_flag)) {
          ++n_vac;
          rho_boost_par = 1;
          rho_wane_par = 1;
          if(sampling_time > vaccination_times[x_vac]) {
            if (sampling_time - vaccination_times[x_vac] > 12) { // within a year 
              wane_amount_vac = 0;
            } else { 
              if (x_vac >= 1) {
                if (vaccination_times[x_vac - 1] == 24180) {
                  rho_boost_par = rho_boost;
                  rho_wane_par = rho_wane;
                }
              }
            }
 
            time = sampling_time - vaccination_times[x_vac]; // Time er vaccination
            wane_amount_vac = MAX(0, 1.0 - (wane_vac*rho_wane_par*time)); // Basic waning function
            seniority = MAX(0, 1.0 - tau_prev_vac*(n_inf + n_vac - 1.0)); // Antigenic seniority
            vac_map_index = vaccination_strain_indices_tmp[x_vac]; // Index of this vaccinating strain in antigenic map

            for (int k = 0; k < n_titres; ++k) {
              index = measurement_strain_indices[tmp_titre_index + k]*number_strains + vac_map_index;
              predicted_titres[tmp_titre_index + k] += (seniority) * ((mu*mu_vac * antigenic_map_long_vac[index]) +
                ((mu_short_vac * rho_boost_par) * antigenic_map_short_vac[index]) * wane_amount_vac);
            }
          }
          ++x_vac;
        }
    }
    start_index_in_data = end_index_in_data;
  }
}