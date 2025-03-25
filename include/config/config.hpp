#pragma once

#include <iostream>
#include <string>
#include <unistd.h>

// store all of our command-line configuration parameters
struct config_t
{
    // The name of the benchmark
    std::string name;

    // The name of the path to the input file
    std::string filename;

    // Number of iterations to run
    int iters;

    // Number of threads
    int threads;

    void print() const {
        std::cout << "# name, filename, iters" << std::endl;
        std::cout << name << ", "
                    << filename << ", "
                    << iters  << ", "
                    << threads << std::endl;
    }
};

// Report on how to use the command line to configure this program
void usage() {
    std::cout
        << "Command-Line Options:" << std::endl
        << "  -n <string> : name of the experiment" << std::endl
        << "  -f <int>    : name of the path to the input file" << std::endl
        << "  -k <int>    : the number of iterations per thread" << std::endl
        << "  -t <int>    : the number of threads in the experiment" << std::endl
        << "  -h          : display this message and exit" << std::endl
        << std::endl;
    exit(0);
}

// Parse command line arguments using getopt()
void parseargs(int argc, char **argv, config_t &cfg) {
    // parse the command-line options
    int opt;
    while ((opt = getopt(argc, argv, "n:f:k:t:h")) != -1) {
        switch (opt) {
        case 'n':
            cfg.name = std::string(optarg);
            break;
        case 'f':
            cfg.filename = std::string(optarg);
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