#ifndef __SRC_LIB_PRICING_HPP__
#define __SRC_LIB_PRICING_HPP__

#include "../vector.hpp"
#include "../matrix.hpp"

class Pricing {  
public:
   virtual HighsInt price(const Vector& x, const Vector& gradient) = 0;
   virtual void update_weights(const Vector& aq, const Vector& ep, HighsInt p, HighsInt q) = 0;
};

#endif
