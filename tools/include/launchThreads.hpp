#pragma once

#include <thread>

// GLobal thread variable for parallelization
int num_threads = 8;

#ifdef SEQUENTIAL_ENABLED
constexpr bool SEQUENTIAL = true;
#else
constexpr bool SEQUENTIAL = false;
#endif

template <typename Func>
void launchThreads(size_t row_num, Func&& f) {
    std::vector<std::thread> threads;
    size_t rows_per_thread = (row_num + num_threads - 1) / num_threads;

    for (int t = 0; t < num_threads; ++t) {
        size_t start = t * rows_per_thread;
        size_t end = std::min(row_num, start + rows_per_thread);
        threads.emplace_back(std::forward<Func>(f), start, end);
    }

    for (auto &thread : threads) {
        thread.join();
    }
}

template <typename Func>
void launchThreadsWithID(size_t row_num, Func&& f) {
    std::vector<std::thread> threads;
    size_t rows_per_thread = (row_num + num_threads - 1) / num_threads;

    for (int t = 0; t < num_threads; ++t) {
        size_t start = t * rows_per_thread;
        size_t end = std::min(row_num, start + rows_per_thread);
        threads.emplace_back(std::forward<Func>(f), start, end, t);
    }

    for (auto &thread : threads) {
        thread.join();
    }
}