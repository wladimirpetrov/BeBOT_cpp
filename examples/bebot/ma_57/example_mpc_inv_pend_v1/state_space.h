#ifndef STATE_SPACE_H
#define STATE_SPACE_H

#include <vector>

std::vector<std::vector<double>> get_A_matrix(double speed);
std::vector<std::vector<double>> get_B_matrix(double speed);

#endif // STATE_SPACE_H
