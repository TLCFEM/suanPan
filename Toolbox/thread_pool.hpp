#pragma once

#include <atomic>
#include <chrono>
#include <condition_variable>
#include <exception>
#include <functional>
#include <future>
#include <iostream>
#include <memory>
#include <mutex>
#include <queue>
#include <thread>
#include <type_traits>
#include <utility>
#include <vector>

/**
 * @file thread_pool.hpp
 * @author Barak Shoshany (baraksh@gmail.com) (http://baraksh.com)
 * @version 3.0.0
 * @date 2022-05-30
 * @copyright Copyright (c) 2022 Barak Shoshany. Licensed under the MIT license. If you use this library in software of any kind, please provide a link to the GitHub repository https://github.com/bshoshany/thread-pool in the source code and documentation. If you use this library in published research, please cite it as follows: Barak Shoshany, "A C++17 Thread Pool for High-Performance Scientific Computing", doi:10.5281/zenodo.4742687, arXiv:2105.00613 (May 2021)
 */

class thread_pool {
    using concurrency_t = std::invoke_result_t<decltype(std::thread::hardware_concurrency)>;

    std::atomic<bool> running = false;

    std::condition_variable task_available_cv = {};

    std::condition_variable task_done_cv = {};

    std::queue<std::function<void()>> tasks = {};

    std::atomic<size_t> tasks_total = 0;

    mutable std::mutex tasks_mutex = {};

    concurrency_t thread_count = 0;

    std::unique_ptr<std::thread[]> threads = nullptr;

    std::atomic<bool> waiting = false;

    void create_threads() {
        running = true;
        for(concurrency_t i = 0; i < thread_count; ++i) threads[i] = std::thread(&thread_pool::worker, this);
    }

    void destroy_threads() {
        running = false;
        task_available_cv.notify_all();
        for(concurrency_t i = 0; i < thread_count; ++i) threads[i].join();
    }

    void worker() {
        while(running) {
            std::function<void()> task;
            std::unique_lock<std::mutex> tasks_lock(tasks_mutex);
            task_available_cv.wait(tasks_lock, [&] { return !tasks.empty() || !running; });
            if(running && !paused) {
                task = std::move(tasks.front());
                tasks.pop();
                tasks_lock.unlock();
                task();
                --tasks_total;
                if(waiting) task_done_cv.notify_one();
            }
        }
    }

public:
    explicit thread_pool(const concurrency_t thread_count_ = std::thread::hardware_concurrency())
        : thread_count(thread_count_ > 1 ? thread_count_ : 1)
        , threads(std::make_unique<std::thread[]>(thread_count)) { create_threads(); }

    ~thread_pool() {
        wait_for_tasks();
        destroy_threads();
    }

    size_t get_tasks_queued() const {
        const std::scoped_lock tasks_lock(tasks_mutex);
        return tasks.size();
    }

    size_t get_tasks_running() const {
        const std::scoped_lock tasks_lock(tasks_mutex);
        return tasks_total - tasks.size();
    }

    size_t get_tasks_total() const { return tasks_total; }

    concurrency_t get_thread_count() const { return thread_count; }

    template<typename F, typename... A> void push_task(const F& task, const A&... args) {
        {
            const std::scoped_lock tasks_lock(tasks_mutex);
            if constexpr(sizeof...(args) == 0) tasks.push(std::function<void()>(task));
            else tasks.push(std::function<void()>([task, args...] { task(args...); }));
        }
        ++tasks_total;
        task_available_cv.notify_one();
    }

    void reset(const concurrency_t thread_count_ = std::thread::hardware_concurrency()) {
        const bool was_paused = paused;
        paused = true;
        wait_for_tasks();
        destroy_threads();
        thread_count = thread_count_ > 1 ? thread_count_ : 1;
        threads = std::make_unique<std::thread[]>(thread_count);
        paused = was_paused;
        create_threads();
    }

    template<typename F, typename... A, typename R = std::invoke_result_t<std::decay_t<F>, std::decay_t<A>...>> std::future<R> submit(const F& task, const A&... args) {
        std::shared_ptr<std::promise<R>> task_promise = std::make_shared<std::promise<R>>();
        push_task(
            [task, args..., task_promise] {
                try {
                    if constexpr(std::is_void_v<R>) {
                        task(args...);
                        task_promise->set_value();
                    }
                    else { task_promise->set_value(task(args...)); }
                }
                catch(...) {
                    try { task_promise->set_exception(std::current_exception()); }
                    catch(...) { }
                }
            });
        return task_promise->get_future();
    }

    void wait_for_tasks() {
        waiting = true;
        std::unique_lock<std::mutex> tasks_lock(tasks_mutex);
        task_done_cv.wait(tasks_lock, [this] { return (tasks_total == (paused ? tasks.size() : 0)); });
        waiting = false;
    }

    std::atomic<bool> paused = false;
};
