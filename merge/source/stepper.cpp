#include "stepper.hpp"

Stepper::Stepper(uint nstages, uint order, double dt)
  : nstages(nstages), irk(0), drk(nstages,0), _t(0.), _dt(dt)
{
  //Allocate the time stepping arrays
  ark.reserve(nstages);
  brk.reserve(nstages);
  crk.reserve(nstages);

  for (uint step = 1; step <= nstages; ++step) {
    ark.emplace_back(std::vector<double>(step,0));
    brk.emplace_back(std::vector<double>(step,0));
    crk.emplace_back(std::vector<double>(step,0));
  }

  //The forward Euler method
  if ((nstages == 1) && (order == 1)) {

    ark[0][0] = 1.;
    brk[0][0] = 1.;

  //SSP(s,2) schemes
  } else if ((nstages == 2) && (order == 2)) {

    ark[0][0] = 1;
    ark[1][0] = 0.5;
    ark[1][1] = 0.5;

    brk[0][0] = 1.;
    brk[1][1] = 0.5;

  //SSP(3,3) scheme
  } else if ((nstages == 3) && (order == 3)) {

    ark[0][0] = 1.;
    ark[1][0] = 3./4.;
    ark[1][1] = 1./4.;
    ark[2][0] = 1./3.;
    ark[2][2] = 2./3.;

    brk[0][0] = 1.;
    brk[1][1] = 1./4.;
    brk[2][2] = 2./3.;

  //SP(4,3) scheme
  } else if ((nstages == 4) && (order == 3)) {

    ark[0][0] = 1.;
    ark[1][1] = 1.;
    ark[2][0] = 2./3.;
    ark[2][2] = 1./3.;
    ark[3][3] = 1.;

    brk[0][0] = 1./2.;
    brk[1][1] = 1./2.;
    brk[2][2] = 1./6.;
    brk[3][3] = 1./2.;

    //SSP(5,3) scheme
  } else if ((nstages == 5) && (order == 3)) {

    ark[0][0] = 1.;
    ark[1][1] = 1.;
    ark[2][0] = 0.355909775063327;
    ark[2][2] = 0.644090224936674;
    ark[3][0] = 0.367933791638137;
    ark[3][3] = 0.632066208361863;
    ark[4][2] = 0.237593836598569;
    ark[4][4] = 0.762406163401431;

    brk[0][0] = 0.377268915331368;
    brk[1][1] = 0.377268915331368;
    brk[2][2] = 0.242995220537396;
    brk[3][3] = 0.238458932846290;
    brk[4][4] = 0.287632146308408;

    //SSP(6,3) scheme
  } else if ((nstages == 6) && (order == 3)) {

    ark[0][0] = 1.;
    ark[1][1] = 1.;
    ark[2][2] = 1.;
    ark[3][0] = 0.476769811285196;
    ark[3][1] = 0.098511733286064;
    ark[3][3] = 0.424718455428740;
    ark[4][4] = 1.;
    ark[5][2] = 0.155221702560091;
    ark[5][5] = 0.844778297439909;

    brk[0][0] = 0.284220721334261;
    brk[1][1] = 0.284220721334261;
    brk[2][2] = 0.284220721334261;
    brk[3][3] = 0.120713785765930;
    brk[4][4] = 0.284220721334261;
    brk[5][5] = 0.240103497065900;

    //SSP(7,3) scheme
  } else if ((nstages == 7) && (order == 3)) {

    ark[0][0] = 1.;
    ark[1][1] = 1.;
    ark[2][2] = 1.;
    ark[3][0] = 0.184962588071072;
    ark[3][3] = 0.815037411928928;
    ark[4][0] = 0.180718656570380;
    ark[4][1] = 0.314831034403793;
    ark[4][4] = 0.504450309025826;
    ark[5][5] = 1.;
    ark[6][3] = 0.120199000000000;
    ark[6][6] = 0.879801000000000;

    brk[0][0] = 0.233213863663009;
    brk[1][1] = 0.233213863663009;
    brk[2][2] = 0.233213863663009;
    brk[3][3] = 0.190078023865845;
    brk[4][4] = 0.117644805593912;
    brk[5][5] = 0.233213863663009;
    brk[6][6] = 0.205181790464579;

    //SSP(8,3) scheme
  } else if ((nstages == 8) && (order == 3)) {

    ark[0][0] = 1.;
    ark[1][1] = 1.;
    ark[2][2] = 1.;
    ark[3][3] = 1.;
    ark[4][0] = 0.421366967085359;
    ark[4][1] = 0.005949401107575;
    ark[4][4] = 0.572683631807067;
    ark[5][1] = 0.004254010666365;
    ark[5][5] = 0.995745989333635;
    ark[6][2] = 0.104380143093325;
    ark[6][3] = 0.243265240906726;
    ark[6][6] = 0.652354615999950;
    ark[7][7] = 1.;

    brk[0][0] = 0.195804015330143;
    brk[1][1] = 0.195804015330143;
    brk[2][2] = 0.195804015330143;
    brk[3][3] = 0.195804015330143;
    brk[4][4] = 0.112133754621673;
    brk[5][5] = 0.194971062960412;
    brk[6][6] = 0.127733653231944;
    brk[7][7] = 0.195804015330143;

    //SSP(5,4) scheme
  } else if ((nstages == 5) && (order == 4)) {

    ark[0][0] = 1.;
    ark[1][0] = 0.44437049406734;
    ark[1][1] = 0.55562950593266;
    ark[2][0] = 0.62010185138540;
    ark[2][2] = 0.37989814861460;
    ark[3][0] = 0.17807995410773;
    ark[3][3] = 0.82192004589227;
    ark[4][0] = 0.00683325884039;
    ark[4][2] = 0.51723167208978;
    ark[4][3] = 0.12759831133288;
    ark[4][4] = 0.34833675773694;

    brk[0][0] = 0.39175222700392;
    brk[1][1] = 0.36841059262959;
    brk[2][2] = 0.25189177424738;
    brk[3][3] = 0.54497475021237;
    brk[4][3] = 0.08460416338212;
    brk[4][4] = 0.22600748319395;

    //SSP(6,4) scheme
  } else if ((nstages == 6) && (order == 4)) {

    ark[0][0] = 1.00000000000000;
    ark[1][0] = 0.30948026455053;
    ark[1][1] = 0.69051973544947;
    ark[2][0] = 0.54205244285557;
    ark[2][2] = 0.45794755714443;
    ark[3][0] = 0.35984960863377;
    ark[3][3] = 0.64015039136623;
    ark[4][4] = 1.00000000000000;
    ark[5][0] = 0.05776282890116;
    ark[5][2] = 0.44216432622405;
    ark[5][4] = 0.10115567086469;
    ark[5][5] = 0.39891717401009;

    brk[0][0] = 0.39270746575722;
    brk[1][1] = 0.30154043149172;
    brk[2][2] = 0.19997937335132;
    brk[3][3] = 0.27954483459696;
    brk[4][4] = 0.43668618869443;
    brk[5][2] = 0.09150931531680;
    brk[5][4] = 0.04417328437472;
    brk[5][5] = 0.14911300530736;

    //SSP(7,4) scheme
  } else if ((nstages == 7) && (order == 4)) {

    ark[0][0] = 1.;
    ark[1][0] = 0.20161507213829;
    ark[1][1] = 0.79838492786171;
    ark[2][0] = 0.19469598207921;
    ark[2][2] = 0.80530401792079;
    ark[3][0] = 0.58143386885601;
    ark[3][3] = 0.41856613114399;
    ark[4][0] = 0.01934367892154;
    ark[4][4] = 0.98065632107846;
    ark[5][5] = 1.;
    ark[6][0] = 0.06006304558847;
    ark[6][2] = 0.30152730794242;
    ark[6][3] = 0.10518998496676;
    ark[6][4] = 0.01483791154585;
    ark[6][6] = 0.51838174995650;

    brk[0][0] = 0.30111872706068;
    brk[1][1] = 0.24040865318216;
    brk[2][2] = 0.24249212077315;
    brk[3][3] = 0.12603810060080;
    brk[4][4] = 0.29529398308716;
    brk[5][5] = 0.30111872706068;
    brk[6][2] = 0.09079551914158;
    brk[6][3] = 0.02888359354880;
    brk[6][6] = 0.15609445267839;

    //SSP(8,4) scheme
  } else if ((nstages == 8) && (order == 4)) {

    ark[0][0] = 1.;
    ark[1][0] = 0.10645325745007;
    ark[1][1] = 0.89354674254993;
    ark[2][2] = 1.;
    ark[3][0] = 0.57175518477257;
    ark[3][3] = 0.42824481522743;
    ark[4][0] = 0.19161667219044;
    ark[4][4] = 0.80838332780956;
    ark[5][5] = 1.;
    ark[6][6] = 1.;
    ark[7][0] = 0.02580435327923;
    ark[7][2] = 0.03629901341774;
    ark[7][3] = 0.31859181340256;
    ark[7][4] = 0.05186768980103;
    ark[7][5] = 0.03944076217320;
    ark[7][6] = 0.00511633747411;
    ark[7][7] = 0.52288003045213;

    brk[0][0] = 0.24120020561311;
    brk[1][1] = 0.21552365802797;
    brk[2][2] = 0.24120020561311;
    brk[3][3] = 0.10329273748560;
    brk[4][4] = 0.19498222488188;
    brk[5][5] = 0.24120020561311;
    brk[6][6] = 0.24120020561311;
    brk[7][2] = 0.00875532949991;
    brk[7][3] = 0.06195575835101;
    brk[7][5] = 0.00951311994571;
    brk[7][7] = 0.12611877085604;

  } else {
    throw std::logic_error("Error invalid Runge-Kutta method entered");
  }

  //Compute the time dependent parameters

  for ( uint i = 0; i < nstages; ++i ) {
    for ( uint k = 0; k < i; ++k ) {
      double casum = 0.;
      for ( uint l = k+1; l < i; ++l ) {
        casum += crk.at(l).at(k)*ark.at(i).at(l);
      }
      crk.at(i).at(k) = brk.at(i).at(k) + casum;
    }
  }

  for( uint k = 1; k < nstages; ++k ) {
    for ( uint l = 0; l < k; ++l ) {
      drk[k] = drk[k] + crk[k-1][l];
    }
  }

  //Compute the maximum beta over alpha ratio at each stage

  /*for ( uint irk = 0; irk < nstages; ++irk ) {
    double max_boa = 0.;
    for ( uint i = 0; i <= irk; ++irk ) {
      if (ark(irk,i) != 0) {
        if (max_boa < brk[irk][i]/ark[irk][i]) { max_boa = brk[irk][i]/ark[irk][i]; }
      }
    }
    max_boa_dt[irk] = max_boa_dt*dt;
    }*/
}