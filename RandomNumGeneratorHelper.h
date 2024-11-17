#include <random>

namespace RandomGenerator{
    int generateRandomNumber(const unsigned int min, const unsigned int max) {
        static std::random_device rd;
        static std::mt19937 gen(rd());

        std::uniform_int_distribution<int> dis(min, max);

        return dis(gen);
    }

    int generateRandomNumber(){
        return generateRandomNumber(0, 100);
    }
}
