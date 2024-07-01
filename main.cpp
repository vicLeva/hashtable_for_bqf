#include <iostream>
#include <fstream>
#include <unordered_map>
#include <string>
#include <sstream>
#include <vector>
#include <ankerl/unordered_dense.h>
#include "additionnal_methods.hpp"


// Function to read k-mers from a file and store them in a hash table
ankerl::unordered_dense::map<u_int64_t, uint8_t> build_index_from_file(const std::string &filename) {
    ankerl::unordered_dense::map<u_int64_t, uint8_t> kmer_table;
    std::ifstream infile(filename);
    std::string line;
    
    while (std::getline(infile, line)) {
        std::istringstream iss(line);
        std::string kmer;
        int count;
        
        if (!(iss >> kmer >> count)) {
            continue; // Skip malformed lines
        }

        if (count > 255) {
            count = 255; // Cap the count at 255 to fit in 8 bits
        }
        
        kmer_table[encode(kmer)] = static_cast<uint8_t>(count);
    }
    
    return kmer_table;
}

// Function to query the hash table with a sequence and get the counts of each k-mer
std::vector<uint8_t> query_kmers(const ankerl::unordered_dense::map<u_int64_t, uint8_t> &kmer_table, const std::string &sequence, int k) {
    std::vector<uint8_t> kmer_counts;
    
    for (size_t i = 0; i <= sequence.length() - k; ++i) {
        uint64_t kmer = canonical(encode(sequence.substr(i, k)), 2*k);
        
        if (kmer_table.find(kmer) != kmer_table.end()) {
            kmer_counts.push_back(kmer_table.at(kmer));
        } else {
            kmer_counts.push_back(0); // If k-mer is not found, count is 0
        }
    }
    
    return kmer_counts;
}



// Function to write the k-mer index to a binary file
void write_index(const std::string &filename, const ankerl::unordered_dense::map<u_int64_t, uint8_t> &kmer_table) {
    std::ofstream outfile(filename, std::ios::binary);
    if (!outfile) {
        throw std::runtime_error("Cannot open output file");
    }

    for (const auto &entry : kmer_table) {
        outfile.write(reinterpret_cast<const char*>(&entry.first), sizeof(entry.first));
        outfile.write(reinterpret_cast<const char*>(&entry.second), sizeof(entry.second));
    }
}

// Function to read the k-mer index from a binary file
ankerl::unordered_dense::map<u_int64_t, uint8_t> read_index(const std::string &filename) {
    ankerl::unordered_dense::map<u_int64_t, uint8_t> kmer_table;
    std::ifstream infile(filename, std::ios::binary);
    if (!infile) {
        throw std::runtime_error("Cannot open input file");
    }

    uint64_t kmer;
    uint8_t count;
    while (infile.read(reinterpret_cast<char*>(&kmer), sizeof(kmer))) {
        infile.read(reinterpret_cast<char*>(&count), sizeof(count));
        kmer_table[kmer] = count;
    }

    return kmer_table;
}

/// Function to handle command line arguments
bool parse_arguments(int argc, char* argv[], std::string &mode, std::vector<std::string> &filenames) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <mode> <files...>" << std::endl;
        std::cerr << "Modes: build <kmer_file> <output_index_file>" << std::endl;
        std::cerr << "       query <input_index_file> <query_file> <output_file>" << std::endl;
        return false;
    }

    mode = argv[1];
    for (int i = 2; i < argc; ++i) {
        filenames.push_back(argv[i]);
    }

    return true;
}

int main(int argc, char* argv[]) {
    std::string mode;
    std::vector<std::string> filenames;
    
    if (!parse_arguments(argc, argv, mode, filenames)) {
        return 1;
    }

    if (mode == "build") {
        if (filenames.size() != 2) {
            std::cerr << "Usage: " << argv[0] << " build <kmer_file> <output_index_file>" << std::endl;
            return 1;
        }

        try {
            ankerl::unordered_dense::map<u_int64_t, uint8_t> kmer_table = build_index_from_file(filenames[0]);
            write_index(filenames[1], kmer_table);
        } catch (const std::exception &e) {
            std::cerr << "Error: " << e.what() << std::endl;
            return 1;
        }

    } else if (mode == "query") {
        if (filenames.size() != 3) {
            std::cerr << "Usage: " << argv[0] << " query <input_index_file> <query_file> <output_file>" << std::endl;
            return 1;
        }

        try {
            ankerl::unordered_dense::map<u_int64_t, uint8_t> kmer_table = read_index(filenames[0]);
            std::ifstream query_file(filenames[1]);
            if (!query_file.is_open()) {
                throw std::runtime_error("Could not open query file: " + filenames[1]);
            }
            
            std::ofstream output_file(filenames[2]);
            if (!output_file.is_open()) {
                throw std::runtime_error("Could not open output file: " + filenames[2]);
            }
            
            std::string line;
            int k = 31; // k-mer length, can be adjusted if needed

            while (std::getline(query_file, line)) {
                if (line.empty() || line[0] == '>') {
                    continue; // Skip header lines and empty lines
                }

                if (line.length() < k) {
                    std::cerr << "Sequence length is less than k: " << line << std::endl;
                    continue;
                }
                
                std::vector<uint8_t> kmer_counts = query_kmers(kmer_table, line, k);
                
                for (uint8_t count : kmer_counts) {
                    output_file << static_cast<int>(count) << " ";
                }
                output_file << std::endl;
            }

        } catch (const std::exception &e) {
            std::cerr << "Error: " << e.what() << std::endl;
            return 1;
        }
    } else {
        std::cerr << "Unknown mode: " << mode << std::endl;
        return 1;
    }

    return 0;
}
