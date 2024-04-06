#include <iostream>
#include <vector>

class MyClass {
private:
    struct MyData {
        int id;
        std::string name;
    };

    std::vector<MyData> dataVector;

    // Private function to generate and return data
    std::vector<MyData> generateData() {
        std::vector<MyData> newData;

        // Generate sample data (can be replaced with actual data generation logic)
        newData.push_back({1, "Data1"});
        newData.push_back({2, "Data2"});
        newData.push_back({3, "Data3"});

        return newData;
    }

public:
    MyClass() {
        // Initialize data using the private function
        dataVector = generateData();
    }

    // Public function to access data
    void printData() {
        for (const auto& item : dataVector) {
            std::cout << "ID: " << item.id << ", Name: " << item.name << std::endl;
        }
    }
};

int main() {
    MyClass myObject;
    myObject.printData();

    return 0;
}

