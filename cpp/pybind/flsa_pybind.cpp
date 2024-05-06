#define PYBIND11_DETAILED_ERROR_MESSAGES 1
#include <vector>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include "flsa.hpp"

namespace py = pybind11;

auto run_flsa(py::array_t<double> y, double lam) {
    // convert NumPy 1D array to C++ vector
    const auto &buf = y.request();
    const auto &shape = buf.shape;
    std::vector<double> y_vec(shape[0]);
    for (int i = 0; i < shape[0]; i++) {
        y_vec[i] = *y.data(i);
    }

    // call the C++ function
    int n = y_vec.size();
    std::vector<int> c(n - 1, 0);
    std::vector<double> beta(n);

    flsa(n, y_vec, lam, c, beta);

    // convert the result back to a NumPy array
    const int size = beta.size();
    py::array_t<double> result(size);
    py::buffer_info numpy_buf = result.request();
    double *ptr = static_cast<double *>(numpy_buf.ptr);
    std::memcpy(ptr, beta.data(), size * sizeof(double));

    return result;
}

PYBIND11_MODULE(tv_denoiser, m) {
    m.def("run_flsa", &run_flsa, "Run FLSA", py::arg("y"), py::arg("lam"));
}
