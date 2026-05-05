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
  if (coeffsin.empty()) {
    this->coeffs = {0.0f};
    return;
  }
  // removing leading 0s
  while (coeffsin.size() > 1 && coeffsin[0] == 0.0f) {
    coeffsin.erase(coeffsin.begin());
  }
  this->coeffs = coeffsin;
}

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

  // TODO: to be implemented

  return roots;
}


std::vector<float> bisection(const Polynomial& p, float tolerance) {
  std::vector<float> roots = {};
  const float cb = p.cauchy_bounds();
  int degree = p.size() - 1;
  
  std::cout << "cauchy bound is: " << cb << std::endl;

  // We will scan the interval [-cb, cb] in small steps to find sign changes
  int num_steps = 1000;
  float step_size = (2.0f * cb) / num_steps;
  
  float current_x = -cb;
  float current_y = p.eval(current_x);
  
  for (int i = 0; i < num_steps; ++i) {
    if (roots.size() == degree) break; // We found all possible real roots

    float next_x = current_x + step_size;
    float next_y = p.eval(next_x);
    
    // Check if a root is exactly at current_x
    if (std::abs(current_y) < tolerance) {
      if (roots.empty() || std::abs(roots.back() - current_x) > tolerance) {
        roots.push_back(current_x);
      }
    } 
    // Check for a sign change indicating a root between current_x and next_x
    else if (std::signbit(current_y) != std::signbit(next_y)) {
      float left = current_x;
      float right = next_x;
      float f_left = current_y;
      
      // Perform bisection in this sub-interval
      while (std::abs(right - left) > tolerance) {
        float mpt = (left + right) / 2.0f;
        float f_mpt = p.eval(mpt);

        if (std::abs(f_mpt) < tolerance) {
          left = mpt;
          right = mpt;
          break;
        }

        if (std::signbit(f_left) != std::signbit(f_mpt)) {
          right = mpt;
        } else {
          left = mpt;
          f_left = f_mpt;
        }
      }
      roots.push_back((left + right) / 2.0f);
    }
    
    current_x = next_x;
    current_y = next_y;
  }
  
  // Check the final endpoint
  if (roots.size() < degree && std::abs(current_y) < tolerance) {
    if (roots.empty() || std::abs(roots.back() - current_x) > tolerance) {
      roots.push_back(current_x);
    }
  }

  if (roots.empty()) {
    std::cout << "a double root situation (y=x^2 kinda thing) , or, no-root situation" << std::endl;
  }

  return roots;
}


float Polynomial::eval(float x) const {
  // evaluates the polynomial at a given point using Horner's method
  if (this->size() == 0) return 0.0f;
  float res = *this->coeff_at(0);
  for (int i = 1; i < this->size(); i++) {
    res = res * x + *this->coeff_at(i);
  }
  return res;
}
