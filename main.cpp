#include <time.h>
#include "LR_lattice.h"
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <sstream>

int main(int argc, char *argv[])
{
    std::cout << "Entering main" << std::endl;
    Init_system();
    SolveEquations(500);
    return 0;
}
