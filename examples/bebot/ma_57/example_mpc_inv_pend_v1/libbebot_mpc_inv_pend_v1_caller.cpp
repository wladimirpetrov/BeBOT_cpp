#include <iostream>
#include <vector>
#include <dlfcn.h>

typedef void (*SolveProblemFunc)(int, double, double, double, double, double, std::vector<double>&);

int main() {
    const char* lib_path = "./libbebot_mpc_inv_pend_v1.so";
    void* handle = dlopen(lib_path, RTLD_LAZY);
    if (!handle) {
        std::cerr << "Failed to load library: " << dlerror() << std::endl;
        return -1;
    }

    SolveProblemFunc solve_problem = (SolveProblemFunc)dlsym(handle, "solve_point_set_problem");
    if (!solve_problem) {
        std::cerr << "Failed to load symbol: " << dlerror() << std::endl;
        dlclose(handle);
        return -1;
    }

    int N = 10;
    double tf = 5.0, theta0 = 1.0, thetadot0 = 0.0, thetaf = 0.0, thetadotf = 0.0;

    std::vector<double> solution;
    solve_problem(N, tf, theta0, thetadot0, thetaf, thetadotf, solution);

    std::cout << "Optimal solution received from .so:" << std::endl;
    for (size_t i = 0; i < solution.size(); ++i) {
        std::cout << "x[" << i << "] = " << solution[i] << std::endl;
    }

    dlclose(handle);
    return 0;
}
