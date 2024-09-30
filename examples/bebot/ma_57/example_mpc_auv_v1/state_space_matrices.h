#ifndef STATE_SPACE_MATRICES_H
#define STATE_SPACE_MATRICES_H

#include <vector>

// Function to get state space matrices A and B based on the speed
void getStateSpaceMatrices(double speed, std::vector<double>& A, std::vector<double>& B);

#endif // STATE_SPACE_MATRICES_H

