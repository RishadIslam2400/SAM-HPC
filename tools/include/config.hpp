#pragma once

#include <iostream>
#include <string>
#include <unistd.h>

// store all our command-line arguments
struct config_t {
    std::string name;        // name of the program
    std::string source_file; // path to the source matrix
    std::string target_file; // path to the taget matrix
    int iters;               // number of iteration for the tests
    int threads;             // specify number of threads

    config_t () {
        name = "default name";
        source_file = "default path";
        target_file = "default path";
        iters = 10;
        threads = 1; // single thread execution by default
    }

    config_t (std::string& progName, std::string& sourceFilePath, std::string& targetFilePath, int numIters, int numThreads) {
        name = progName;
        source_file = sourceFilePath;
        target_file = targetFilePath;
        iters = numIters;
        threads = numThreads;
    }

    void print() {
        std::cout << "# name, sourcefile, targetfile, iters, threads:\n";
        std::cout << name << ", " << source_file << ", " << target_file << ", " << iters << ", " << threads << std::endl;
    }
};

// report on how to use the command line to configure this program
void usage() {
    std::cout << "Command line options:\n"
              << "-n <string> : name of the experiment\n"
              << "-x <string> : path to the source matrix file\n"
              << "-y <string> : path to the target matrix file\n" 
              << "-k <int> : number of test iterations\n"
              << "-t : number of threads for the execution\n"
              << "-h : display help message"
              << std::endl;
}

// parse command line arguments using get-opt()
void parseargs(int argc, char** argv, config_t& cfg) {
    int opt;
    while ((opt = getopt(argc, argv, "n:x:y:k:t:h")) != -1) {
        switch (opt) {
        case 'n':
            cfg.name = std::string(optarg);
            break;
        case 'x':
            cfg.source_file = std::string(optarg);
            break;
        case 'y':
            cfg.target_file = std::string(optarg);
            break;
        case 'k':
            cfg.iters = atoi(optarg);
            break;
        case 't':
            cfg.threads = atoi(optarg);
            break;
        case 'h':
            usage();
            break;
        }
    }
}