/**
 * @class thread_pool
 * @brief A fixed-size thread pool for asynchronous task execution.
 *
 * Abridged from BS::thread_pool and limited to the features used by suanPan.
 *
 * @author tlc
 * @date 15/05/2022
 * @version 0.1.0
 * @file thread_pool.hpp
 * @addtogroup Utility
 * @{
 */
#ifndef THREAD_POOL_HPP
#define THREAD_POOL_HPP

#include <condition_variable>
#include <future>
#include <queue>
#include <suanPan.h>

class [[nodiscard]] thread_pool {
    std::mutex tasks_mutex;

    std::condition_variable task_available_cv;

    std::condition_variable tasks_done_cv;

    std::queue<std::function<void()>> tasks;

    // Number of tasks currently executing, guarded by tasks_mutex.
    std::size_t tasks_running = 0;

    // Number of threads blocked in wait_for_tasks(), guarded by tasks_mutex.
    std::size_t waiters = 0;

    bool running = false;

    std::size_t thread_count = 0;

    std::vector<std::thread> threads;

    [[nodiscard]] static std::size_t determine_thread_count(const std::size_t count) {
        if(count > 0) return count;
        const auto available = static_cast<std::size_t>(std::thread::hardware_concurrency());
        return 0 == available ? 1 : available;
    }

    void create_threads() {
        {
            const std::scoped_lock lock(tasks_mutex);
            running = true;
        }

        try {
            threads.reserve(thread_count);
            for(std::size_t I = 0; I < thread_count; ++I) threads.emplace_back(&thread_pool::worker, this);
        }
        catch(...) {
            {
                const std::scoped_lock lock(tasks_mutex);
                running = false;
            }
            task_available_cv.notify_all();
            for(auto& thread : threads) thread.join();
            throw;
        }
    }

    void destroy_threads() {
        {
            const std::scoped_lock lock(tasks_mutex);
            running = false;
        }
        task_available_cv.notify_all();
        for(auto& thread : threads) thread.join();
    }

    void worker() {
        while(true) {
            std::function<void()> task;
            {
                std::unique_lock lock(tasks_mutex);
                task_available_cv.wait(lock, [this] { return !running || !tasks.empty(); });
                if(!running) return;
                task = std::move(tasks.front());
                tasks.pop();
                ++tasks_running;
            }

            try {
                task();
            }
            catch(const std::exception& error) {
                try {
                    suanpan_error("Unhandled exception in thread pool task: {}\n", error.what());
                }
                catch(...) {
                }
            }
            catch(...) {
                try {
                    suanpan_error("Unhandled non-standard exception in thread pool task.\n");
                }
                catch(...) {
                }
            }

            {
                const std::scoped_lock lock(tasks_mutex);
                --tasks_running;
                if(waiters > 0 && 0 == tasks_running && tasks.empty()) tasks_done_cv.notify_all();
            }
        }
    }

    template<typename F> void enqueue(F&& task) {
        {
            const std::scoped_lock lock(tasks_mutex);
            tasks.emplace(std::forward<F>(task));
        }
        task_available_cv.notify_one();
    }

public:
    explicit thread_pool(const std::size_t count = 0)
        : thread_count(determine_thread_count(count)) { create_threads(); }

    thread_pool(const thread_pool&) = delete;
    thread_pool(thread_pool&&) = delete;
    thread_pool& operator=(const thread_pool&) = delete;
    thread_pool& operator=(thread_pool&&) = delete;

    ~thread_pool() {
        wait_for_tasks();
        destroy_threads();
    }

    [[nodiscard]] std::size_t get_thread_count() const { return thread_count; }

    template<typename F, typename... A> void push_task(F&& task, A&&... args) {
        // Skip the indirection for the common argument-less task.
        if constexpr(0 == sizeof...(A)) enqueue(std::forward<F>(task));
        else enqueue([task = std::forward<F>(task), ... args = std::forward<A>(args)]() mutable { std::invoke(std::move(task), std::move(args)...); });
    }

    template<typename F, typename... A, typename R = std::invoke_result_t<std::decay_t<F>, std::decay_t<A>...>> [[nodiscard]] std::future<R> submit(F&& task, A&&... args) {
        auto packaged = std::make_shared<std::packaged_task<R()>>(
            [task = std::forward<F>(task), ... args = std::forward<A>(args)]() mutable -> R { return std::invoke(std::move(task), std::move(args)...); }
        );
        auto future = packaged->get_future();
        enqueue([packaged] { (*packaged)(); });
        return future;
    }

    template<typename T1, typename T2, typename F, typename T = std::common_type_t<T1, T2>> void push_loop(const T1 first_index, const T2 index_after_last, F&& loop, const std::size_t num_blocks = 0) {
        auto first = static_cast<T>(first_index);
        auto last = static_cast<T>(index_after_last);
        if(last < first) std::swap(first, last);

        const auto total_size = static_cast<std::size_t>(last - first);
        if(0 == total_size) return;

        auto blocks = 0 == num_blocks ? thread_count : num_blocks;
        auto block_size = total_size / blocks;
        if(0 == block_size) {
            block_size = 1;
            blocks = total_size;
        }

        const auto block_loop = std::forward<F>(loop);
        for(std::size_t I = 0; I < blocks; ++I) push_task(block_loop, static_cast<T>(I * block_size) + first, I == blocks - 1 ? last : static_cast<T>((I + 1) * block_size) + first);
    }

    template<typename T, typename F> void push_loop(const T index_after_last, F&& loop, const std::size_t num_blocks = 0) { push_loop(T{0}, index_after_last, std::forward<F>(loop), num_blocks); }

    void wait_for_tasks() {
        std::unique_lock lock(tasks_mutex);
        if(0 == tasks_running && tasks.empty()) return;
        ++waiters;
        tasks_done_cv.wait(lock, [this] { return 0 == tasks_running && tasks.empty(); });
        --waiters;
    }
};

#endif

//! @}
