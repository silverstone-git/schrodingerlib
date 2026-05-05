#ifndef POLYNOMIAL_H

#define POLYNOMIAL_H

#include <optional>
#include <vector>
class Polynomial {
  public:
    explicit Polynomial(std::vector<float> coeff_vectors);
    std::optional<float> coeff_at(int i) const;
    int size() const;
    float eval(float x) const;
    float cauchy_bounds() const;
    std::vector<float> roots() const;
    void print() const;

  private:
    std::vector<float> coeffs;
};

std::vector<float> bisection(const Polynomial& p, float tolerance = 0.0001f);

#endif
