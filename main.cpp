#include "tvector.h"

using namespace std;

int main()
{
    TVector vector_1({1, 2, 3});
    TVector vector_2({2, 3, 4});
    // vector_1 = vector_2;
    TVector vector_3 = vector_1 - vector_2;
    // TVector vector_4(vector_2);
    std::cout << "V3 = " << vector_3 << std::endl;
    std::cout << "V1 - V2 = " << vector_1 - vector_2 << std::endl;
    std::cout << "-V1 = " << -vector_1 << std::endl;
    TVector vector_4 = vector_1 * 5;
    std::cout << "V1 = " << vector_1 << std::endl;
    std::cout << "V4* = " << vector_4 << std::endl;
    std::cout << "V4 * V1 = " << vector_4 * vector_1 << std::endl;
    return 0;
}
