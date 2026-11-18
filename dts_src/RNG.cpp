
#include "RNG.h"
// Constructor
RNG::RNG(int seed) : m_seed(seed) {

   m_generator.seed(m_seed);  // Use the seed to initialize the generator
   m_uniform_dist = std::uniform_real_distribution<double>(0.0, 1.0);  // Uniform distribution between 0 and 1
    
    m_distribution = std::normal_distribution<double> (0, 1);

}
int RNG::StreamSeed(int seed, int stream) {
    if (stream == 0) {
        return seed;
    }
    uint64_t z = ((uint64_t)(uint32_t)seed << 32) + (uint64_t)stream * 0x9E3779B97F4A7C15ULL;
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    z ^= z >> 31;
    return (int)(z & 0x7FFFFFFF);
}
void RNG::SetGaussianDistribution(double mean, double stddev) {
    m_distribution = std::normal_distribution<double> (mean, stddev);
    return;
}
// Uniform RNG generator for doubles in range [0, A)
double RNG::UniformRNG(double A) {
        return A * m_uniform_dist(m_generator);  // Use C++11 uniform distribution
}

// Uniform RNG generator for integers in range [0, a)
int RNG::IntRNG(int a) {
        
        std::uniform_int_distribution<int> int_dist(0, a - 1);
        return int_dist(m_generator);  // Use C++11 uniform int distribution
}

// Bernoulli RNG for a binary outcome (true/false)
bool RNG::BinRNG() {
    // Generate a uniform random number in range [0, 1] and return true if greater than 0.5
    return UniformRNG(1.0) > 0.5;
}
double RNG::GaussianRNG() {

    return m_distribution(m_generator);
}
