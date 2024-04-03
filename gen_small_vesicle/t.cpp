#include <iostream>
#include <chrono>
#include <thread>

int main() {
    for (int i = 0; i <= 100; ++i) {
        // Print progress percentage
        std::cout << "Progress: " << i << "%" << std::flush;

        // Simulate some work
        std::this_thread::sleep_for(std::chrono::milliseconds(100));

        // Move the cursor back to the beginning of the line
        std::cout << '\r';

        // Clear the line
        std::cout << "\033[K";
    }

    // Print completion message
    std::cout << "Task completed!" << std::endl;

    return 0;
}
