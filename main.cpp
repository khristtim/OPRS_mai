#include "tvector.h"

using namespace std;

int main()
{
    TVector vector_1({1, 2, 3});
    TVector vector_2({2, 3, 4});
    vector_1 = vector_2;
    // TVector vector_3 = vector_1 - vector_2;
    TVector vector_3(vector_2);
    std::cout << "V1" << vector_3 << std::endl;
    // std::cout << "V2" << vector_1 - vector_2 << std::endl;
    std::cout << "V4" << -vector_1 << std::endl;
    std::cout << "V3" << vector_1 << std::endl;
    return 0;
}
