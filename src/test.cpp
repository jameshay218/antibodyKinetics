#include <Rcpp.h>
using namespace Rcpp;

//' @export
//' @useDynLib antibodyKinetics
//[[Rcpp::export]]
NumericMatrix predict_titres(NumericVector params, NumericVector times){
  double y0 = 0;
  double final_t;
  double t_i;
  double mu, dRF, tp, ts, m;
  double tmp;
  double tmp2;
  NumericVector::iterator t = times.begin();
  NumericVector Y(times.size());
  

  int j = 0;
  double trunc_lower = params[0];
  double first_infection = params[1];
  double max_t = times[times.size()-1];
  int no_infections=params.size()/6;
  NumericMatrix out(times.size(),2);

  for(int i=1;i <= no_infections; ++i){
    if(i==no_infections){
      final_t = max_t;
    }
    else {
      final_t = params[6*(i)+1];
    }

    t_i = params[6*(i-1)+1];
    mu = params[6*(i-1)+2];
    dRF = params[6*(i-1)+3];
    tp = params[6*(i-1)+4];
    ts = params[6*(i-1)+5];
    m = params[6*(i-1)+6];

    while(t != times.end() && *t <= final_t){
      tmp = 0;
      if(*t <= t_i) tmp = 0;
      else if(*t > t_i && *t <= (t_i+tp)) tmp = (mu/tp)*(*t-t_i);
      else if(*t > (tp+t_i) && *t <=(ts + t_i+tp)) tmp = ((-(dRF*mu)/ts)*(*t) + ((mu*dRF)/ts)*(t_i+tp) + mu);
      else tmp = (-m*(*t)+m*(t_i+tp+ts)+(1-dRF)*mu);
      tmp += y0;
      
      Y[j] = tmp;
      ++t;
      ++j;
    }
    
    if((6*i + 1) < params.size()){
      tmp2 = params[6*i + 1];
      if(tmp2 <= t_i) tmp = 0;
      else if(tmp2 > t_i && tmp2 <= (t_i+tp)) tmp = (mu/tp)*(tmp2-t_i);
      else if(tmp2 > (tp+t_i) && tmp2 <=(ts + t_i+tp)) tmp = ((-(dRF*mu)/ts)*(tmp2) + ((mu*dRF)/ts)*(t_i+tp) + mu);
      else tmp = (-m*(tmp2)+m*(t_i+tp+ts)+(1-dRF)*mu);
      y0 = y0 + tmp;
    }
    
    if(y0 < trunc_lower){
      y0 = trunc_lower;
    }
  }
  
  out(_,0) = times;
  out(_,1)=Y;
  return out;
}



//' Converts to unit scale
//'
//' @param x the double to be converted
//' @param min the minimum value on the linear scale
//' @param max the maximum value on the linear scale
//' @return the value converted to a unit scale
//' @export
//[[Rcpp::export]]
double toUnitScale(double x, double min, double max){
  return((x-min)/(max-min));
}

//' Converts to linear scale
//'
//' @param x the double to be converted back to linear scale
//' @param min the minimum value on the linear scale
//' @param max the maximum value on the linear scale
//' @return the value converted to a linear scale
//' @export
//[[Rcpp::export]]
double fromUnitScale(double x, double min, double max){
  return(min + (max-min)*x);
}


//' @export
//[[Rcpp::export]]
double obs_error(int actual, int obs, double S, double EA){
  int MAX_TITRE = 13;
  if(actual == (MAX_TITRE) && obs == (MAX_TITRE)) return(S + EA/2.0 - (1.0/(MAX_TITRE-2.0))*(1.0-S-EA));
  else if(actual==0 && obs==0) return(S + EA/2.0 - (1.0/(MAX_TITRE-2.0))*(1.0-S-EA));
  else if(actual==obs) return(S);
  else if(actual == (obs + 1) || actual==(obs-1)) return(EA/2.0);
  return((1.0/(MAX_TITRE-2.0))*(1.0-S-EA));
}


//' @export
//[[Rcpp::export]]
double posterior(NumericVector y, NumericVector data, NumericVector params){
  double ln = 0;
  for(int i = 0; i < y.length();++i){
    if(y(i) < 0) y(i) = 0;
    if(y(i) >= 13) y(i) = 13;
    ln += log(obs_error(floor(y(i)), floor(data(i)),params(0),params(1)));
  }
  return ln;
}
