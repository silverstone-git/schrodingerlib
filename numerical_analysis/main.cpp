#include "polynomial/Polynomial.h"
#include <iostream>
#include <ostream>
#include <vector>
#include <string>


int main() {

  std::cout << "Enter size of polynomial: ";
  int size;
  std::cin >> size;

  auto MAX_SIZE = 1000;

  if(size < 1 || size > MAX_SIZE) {
    std::cout << "Invalid size!" << std::endl;
    return 1;
  }

  std::string preview_string = "";

  for(int i = 0; i < size; i++) {
    preview_string += "a" + std::to_string(i);
    preview_string += "_x^";
    preview_string += std::to_string(size - i - 1);
    if (i < size - 1) {
      preview_string += " + ";
    }
  }
  std::cout << "Enter Polynomial as coefficients : " << preview_string << " : ";

  std::vector<float> polyn(size);
  for(int i = 0; i < size; i ++) {
    std::cin >> polyn[i];
  }

  Polynomial p = Polynomial(polyn);
  p.print();
  // std::cout << p.eval(1) << std::endl;
  // std::cout << p.eval(2) << std::endl;

  std::vector<float> res = bisection(p);
  for(const auto& root : res) std::cout << std::to_string(root) << ", ";
  std::cout << std::endl;
  return 0;
}
