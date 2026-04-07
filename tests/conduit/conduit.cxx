#include <iostream>
#include "conduit.hpp"
#include "conduit_relay.hpp"

using namespace conduit;

int main(int argc, char** argv) {
    // Check command line arguments
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <input_yaml_file>" << std::endl;
        return 1;
    }
    
    std::string yaml_file = argv[1];
    
    try {
        // Create a Conduit node (data structure)
        Node data;
        
        // Load YAML file into Conduit node
        relay::io::load(yaml_file, "yaml", data);
        
        std::cout << "=== YAML File Contents ===" << std::endl;
        
        // Option 1: Print in YAML format
        std::cout << "\n1. Pretty YAML output:" << std::endl;
        data.print();
        
          
        // Check if a specific path exists
        if (data.has_path("name")) {
            std::string name = data["name"].as_string();
            std::cout << "Name: " << name << std::endl;
        }
        
        // Check for nested values
        if (data.has_path("user/profile/age")) {
            int age = data["user/profile/age"].to_int32();
            std::cout << "Age: " << age << std::endl;
        }
        
        // Handle arrays/lists
        if (data.has_path("scores")) {
            std::cout << "Scores: ";
            Node &scores = data["scores"];
            for (int i = 0; i < scores.dtype().number_of_elements(); i++) {
                std::cout << scores.as_int32_accessor()[i] << " ";
            }
            std::cout << std::endl;
        }
        
    } catch (const conduit::Error &e) {
        std::cerr << "Conduit error: " << e.what() << std::endl;
        return 1;
    } catch (const std::exception &e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}