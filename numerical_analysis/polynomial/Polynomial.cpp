#include "Polynomial.h"
#include <cfloat>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <optional>
#include <ostream>
#include <string>
#include <vector>

Polynomial::Polynomial(std::vector<float> coeffsin) {
  // removing leading 0s
  while (coeffsin[0] == 0.0f && coeffsin.size() > 1) {
    coeffsin.erase(coeffsin.begin());
  }
  this->coeffs = coeffsin;

};

void Polynomial::print() const {

  // std::cout << "{ ";
  // for (const auto& vec_item : this->coeffs) {
  //   std::cout << vec_item << " ";
  // }
  //
  // std::cout << "}" << std::endl;
  
  std::string tbp = "";

  for(int i = 0; i < this->size() - 1; i++) {
    int mirror_index = this->size() - i - 1;
    if(mirror_index == 1) {
      tbp += std::to_string(*this->coeff_at(i)) + "x + ";
    } else {
      tbp += std::to_string(*this->coeff_at(i));
      tbp += "x^";
      tbp += std::to_string(mirror_index);
      tbp += " + ";

    }
  }
  tbp += std::to_string(*this->coeff_at(this->size() - 1));
  std::cout << tbp << std::endl;
}

int Polynomial::size() const {
  return this->coeffs.size();
}

std::optional<float> Polynomial::coeff_at(int i) const {
  if(i > this->size() - 1 || i < 0) {
    return std::nullopt;
  } else {
    return this->coeffs.at(i);
  }
}

float Polynomial::cauchy_bounds() const {
  if(this->coeff_at(0) == 0) {
    throw "Invalid Polynomial passed into bisection: got leading 0";
  }
  float max_ratio = FLT_MIN;
  for(int i = 1; i < this->size(); i++) {
    std::optional<float> cur_coeff = this->coeff_at(i);
    std::optional<float> first_coeff = this->coeff_at(0);
    if(cur_coeff && first_coeff) {
      max_ratio = std::max(max_ratio, std::abs(*cur_coeff / *first_coeff));
    }
  }
  return 1 + max_ratio;
}

auto Polynomial::roots() const -> std::vector<float> {

  std::vector<float> roots = {};

  // edge cases
  if(this->size() == 1) {
    if(this->coeff_at(0) != 0.0f) {
      return roots;
    } else {
      std::cout << "Infinite Roots!" << std::endl;
      return roots;
    }
  }

  return roots;
}


std::vector<float> bisection(const Polynomial& p) {
  const int cauchy_bounds = p.cauchy_bounds();
  int left = -cauchy_bounds;
  int right = cauchy_bounds;

  // while ()
  std::vector<float> roots = {1.0f};
  return roots;
}


float Polynomial::eval(float x) const {
  // evalutes the polynomial at a given point
  int res = 0;
  for(int i = 0; i < this->size() - 1; i ++) {
    res += *coeff_at(i) * pow(x, this->size() - i - 1);
  }
  res += *this->coeff_at(this->size() - 1);
  return res;
}
