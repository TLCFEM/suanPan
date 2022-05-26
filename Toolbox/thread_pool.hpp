// ReSharper disable CommentTypo
#ifndef THREAD_POOL_HPP
#define THREAD_POOL_HPP

/**
 * @file thread_pool.hpp
 * @author Barak Shoshany (baraksh@gmail.com) (http://baraksh.com)
 * @version 2.0.0
 * @date 2021-08-14
 * @copyright Copyright (c) 2021 Barak Shoshany. Licensed under the MIT license. If you use this library in published research, please cite it as follows:
 *  - Barak Shoshany, "A C++17 Thread Pool for High-Performance Scientific Computing", doi:10.5281/zenodo.4742687, arXiv:2105.00613 (May 2021)
 *
 * @brief A C++17 thread pool for high-performance scientific computing.
 * @details A modern C++17-compatible thread pool implementation, built from scratch with high-performance scientific computing in mind. The thread pool is implemented as a single lightweight and self-contained class, and does not have any dependencies other than the C++17 standard library, thus allowing a great degree of portability. In particular, this implementation does not utilize OpenMP or any other high-level multithreading APIs, and thus gives the programmer precise low-level control over the details of the parallelization, which permits more robust optimizations. The thread pool was extensively tested on both AMD and Intel CPUs with up to 40 cores and 80 threads. Other features include automatic generation of futures and easy parallelization of loops. Two helper classes enable synchronizing printing to an output stream by different threads and measuring execution time for benchmarking purposes. Please visit the GitHub repository at https://github.com/bshoshany/thread-pool for documentation and updates, or to submit feature requests and bug reports.
 */

#include <atomic>
#include <chrono>
#include <cstdint>
#include <functional>
#include <future>
#include <memory>
#include <mutex>
#include <queue>
#include <thread>
#include <type_traits>
#include <utility>

class thread_pool {
    using ui32 = std::uint_fast32_t;
    using ui64 = std::uint_fast64_t;

    mutable std::mutex queue_mutex = {};

    std::atomic<bool> running = true;

    std::queue<std::function<void()>> tasks = {};

    ui32 thread_count;

    std::unique_ptr<std::thread[]> threads;

    std::atomic<ui32> tasks_total = 0;

    void create_threads() { for(ui32 i = 0; i < thread_count; i++) threads[i] = std::thread(&thread_pool::worker, this); }

    void destroy_threads() { for(ui32 i = 0; i < thread_count; i++) threads[i].join(); }

    bool pop_task(std::function<void()>& task) {
        const std::scoped_lock lock(queue_mutex);
        if(tasks.empty()) return false;
        task = std::move(tasks.front());
        tasks.pop();
        return true;
    }

    void sleep_or_yield() const {
        if(sleep_duration) std::this_thread::sleep_for(std::chrono::microseconds(sleep_duration));
        else std::this_thread::yield();
    }

    void worker() {
        while(running)
            if(std::function<void()> task; !paused && pop_task(task)) {
                task();
                --tasks_total;
            }
            else sleep_or_yield();
    }

public:
    std::atomic<bool> paused = false;

    ui32 sleep_duration = 1000;

    explicit thread_pool(const ui32 _thread_count = std::thread::hardware_concurrency())
        : thread_count(_thread_count ? _thread_count : std::thread::hardware_concurrency())
        , threads(new std::thread[thread_count]) { create_threads(); }

    thread_pool(const thread_pool&) = delete;
    thread_pool(thread_pool&&) noexcept = delete;
    thread_pool& operator=(const thread_pool&) = delete;
    thread_pool& operator=(thread_pool&&) noexcept = delete;

    ~thread_pool() {
        wait_for_tasks();
        running = false;
        destroy_threads();
    }

    ui64 get_tasks_queued() const {
        const std::scoped_lock lock(queue_mutex);
        return tasks.size();
    }

    ui32 get_tasks_running() const { return tasks_total - static_cast<ui32>(get_tasks_queued()); }

    ui32 get_tasks_total() const { return tasks_total; }

    ui32 get_thread_count() const { return thread_count; }

    template<typename F> void push_task(const F& task) {
        ++tasks_total;
        {
            const std::scoped_lock lock(queue_mutex);
            tasks.push(std::function<void()>(task));
        }
    }

    template<typename F, typename... A> void push_task(const F& task, const A&... args) { push_task([task, args...] { task(args...); }); }

    void reset(const ui32 _thread_count = std::thread::hardware_concurrency()) {
        const bool was_paused = paused;
        paused = true;
        wait_for_tasks();
        running = false;
        destroy_threads();
        thread_count = _thread_count ? _thread_count : std::thread::hardware_concurrency();
        threads.reset(new std::thread[thread_count]);
        paused = was_paused;
        running = true;
        create_threads();
    }

    template<typename F, typename... A, typename = std::enable_if_t<std::is_void_v<std::invoke_result_t<std::decay_t<F>, std::decay_t<A>...>>>> std::future<bool> submit(const F& task, const A&... args) {
        std::shared_ptr<std::promise<bool>> task_promise(new std::promise<bool>);
        std::future<bool> future = task_promise->get_future();
        push_task([task, args..., task_promise] {
            try {
                task(args...);
                task_promise->set_value(true);
            }
            catch(...) {
                try { task_promise->set_exception(std::current_exception()); }
                catch(...) { }
            }
        });
        return future;
    }

    template<typename F, typename... A, typename R = std::invoke_result_t<std::decay_t<F>, std::decay_t<A>...>, typename = std::enable_if_t<!std::is_void_v<R>>> std::future<R> submit(const F& task, const A&... args) {
        std::shared_ptr<std::promise<R>> task_promise(new std::promise<R>);
        std::future<R> future = task_promise->get_future();
        push_task([task, args..., task_promise] {
            try { task_promise->set_value(task(args...)); }
            catch(...) {
                try { task_promise->set_exception(std::current_exception()); }
                catch(...) { }
            }
        });
        return future;
    }

    void wait_for_tasks() const {
        while(true) {
            if(!paused) { if(tasks_total == 0) break; }
            else { if(get_tasks_running() == 0) break; }
            sleep_or_yield();
        }
    }
};

#endif
