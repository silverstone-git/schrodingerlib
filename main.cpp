#include <cmath>
#include <matplot/matplot.h>
#include <vector>

int main() {
    using namespace matplot;
    std::vector<double> x = linspace(0, 2 * pi);
    std::vector<double> y = transform(x, [](auto x) {return sin(x);});
    plot(x, y);

      show();
    return 0;
}
