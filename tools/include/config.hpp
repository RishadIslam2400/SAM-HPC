#pragma once

#include <iostream>
#include <string>
#include <unistd.h>

// store all our command-line arguments
struct config_t {
    std::string name;        // name of the program
    std::string matrix_dir; // path to the source matrix
    std::string rhs_dir; // path to the taget matrix
    int iters;               // number of iteration for the tests
    int threads;             // specify number of threads

    config_t () {
        name = "default name";
        matrix_dir = "default path";
        rhs_dir = "default path";
        iters = 10;
        threads = 1; // single thread execution by default
    }

    config_t (std::string& progName, std::string& matrixDirPath, std::string& rhsDirPath, int numIters, int numThreads) {
        name = progName;
        matrix_dir = matrixDirPath;
        rhs_dir = rhsDirPath;
        iters = numIters;
        threads = numThreads;
    }

    void print() {
        std::cout << "# name, matrix_dir, rhs_dir, iters, threads:\n";
        std::cout << name << ", " << matrix_dir << ", " << rhs_dir << ", " << iters << ", " << threads << std::endl;
    }
};

// report on how to use the command line to configure this program
void usage() {
    std::cout << "Command line options:\n"
              << "-n <string> : name of the experiment\n"
              << "-x <string> : path to the matrix directory\n"
              << "-y <string> : path to the rhs directory\n" 
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
            cfg.matrix_dir = std::string(optarg);
            break;
        case 'y':
            cfg.rhs_dir = std::string(optarg);
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