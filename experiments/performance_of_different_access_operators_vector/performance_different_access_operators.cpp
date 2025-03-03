#include <fmt/format.h>  // For formatted output

#include <chrono>  // For timing
#include <cstdint>
#include <iostream>
#include <vector>

int main()
{
    const uint64_t num_samples = 1000000000;

    const std::vector<uint64_t> values = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};

    // Volatile variable to prevent optimization
    volatile uint64_t temp = 0;

    // Measure time for [] operator
    auto start = std::chrono::high_resolution_clock::now();
    for (uint64_t i = 0; i < num_samples; i++)
    {
        for (uint64_t l = 0; l < values.size(); l++)
        {
            temp += values[l];  // Access elements using []
        }
    }
    auto                                end            = std::chrono::high_resolution_clock::now();
    const std::chrono::duration<double> duration_first = end - start;  // Time for [] operator

    // Measure time for .at() operator
    start = std::chrono::high_resolution_clock::now();
    for (uint64_t i = 0; i < num_samples; i++)
    {
        for (uint64_t l = 0; l < values.size(); l++)
        {
            temp += values.at(l);  // Access elements using .at()
        }
    }
    end                                                 = std::chrono::high_resolution_clock::now();
    const std::chrono::duration<double> duration_second = end - start;  // Time for .at() operator

    // Print the results
    std::cout << fmt::format("Time for [] operator: {:.6f} seconds\n", duration_first.count());
    std::cout << fmt::format("Time for .at() operator: {:.6f} seconds\n", duration_second.count());

    return 0;
}
