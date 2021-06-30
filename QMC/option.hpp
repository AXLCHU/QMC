#ifndef __OPTION_HPP
#define __OPTION_HPP

#include "PayOff.hpp"

class Option {
 public:
  PayOff* pay_off;
  double K;
  double r;
  double T;

  Option(double _K, double _r, double _T, PayOff* _pay_off);

  virtual ~Option();
};

#endif