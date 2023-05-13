#ifndef THREAD_POOL_HPP
#define THREAD_POOL_HPP

#include <atomic>
#include <condition_variable>
#include <exception>
#include <functional>
#include <future>
#include <memory>
#include <mutex>
#include <queue>
#include <thread>
#include <type_traits>
#include <utility>

class [[nodiscard]] thread_pool {
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
        {
            const std::scoped_lock tasks_lock(tasks_mutex);
            task_available_cv.notify_all();
        }
        for(concurrency_t i = 0; i < thread_count; ++i) threads[i].join();
    }

    [[nodiscard]] static concurrency_t determine_thread_count(const concurrency_t thread_count_) { return thread_count_ > 0 ? thread_count_ : std::thread::hardware_concurrency() > 0 ? std::thread::hardware_concurrency() : 1; }

    void worker() {
        while(running) {
            std::unique_lock tasks_lock(tasks_mutex);
            task_available_cv.wait(tasks_lock, [this] { return !tasks.empty() || !running; });
            if(running) {
                std::function<void()> task;
                task = std::move(tasks.front());
                tasks.pop();
                tasks_lock.unlock();
                task();
                tasks_lock.lock();
                --tasks_total;
                if(waiting) task_done_cv.notify_one();
            }
        }
    }

public:
    explicit thread_pool(const concurrency_t thread_count_ = 0)
        : thread_count(determine_thread_count(thread_count_))
        , threads(std::make_unique<std::thread[]>(determine_thread_count(thread_count_))) { create_threads(); }

    ~thread_pool() {
        wait_for_tasks();
        destroy_threads();
    }

    [[nodiscard]] concurrency_t get_thread_count() const { return thread_count; }

    template<typename F, typename T1, typename T2, typename T = std::common_type_t<T1, T2>> void push_loop(T1 first_index_, T2 index_after_last_, F&& loop, size_t num_blocks = 0) {
        T first_index = static_cast<T>(first_index_);
        T index_after_last = static_cast<T>(index_after_last_);
        if(num_blocks == 0) num_blocks = thread_count;
        if(index_after_last < first_index) std::swap(index_after_last, first_index);
        const size_t total_size = static_cast<size_t>(index_after_last - first_index);
        size_t block_size = total_size / num_blocks;
        if(block_size == 0) {
            block_size = 1;
            num_blocks = (total_size > 1) ? total_size : 1;
        }
        if(total_size > 0) { for(size_t i = 0; i < num_blocks; ++i) push_task(std::forward<F>(loop), static_cast<T>(i * block_size) + first_index, (i == num_blocks - 1) ? index_after_last : (static_cast<T>((i + 1) * block_size) + first_index)); }
    }

    template<typename F, typename T> void push_loop(const T index_after_last, F&& loop, const size_t num_blocks = 0) { push_loop(0, index_after_last, std::forward<F>(loop), num_blocks); }

    template<typename F, typename... A> void push_task(F&& task, A&&... args) {
        {
            const std::function<void()> task_function = std::bind(std::forward<F>(task), std::forward<A>(args)...);
            const std::scoped_lock tasks_lock(tasks_mutex);
            tasks.push(task_function);
            ++tasks_total;
        }
        task_available_cv.notify_one();
    }

    template<typename F, typename... A, typename R = std::invoke_result_t<std::decay_t<F>, std::decay_t<A>...>> [[nodiscard]] std::future<R> submit(F&& task, A&&... args) {
        std::function<R()> task_function = std::bind(std::forward<F>(task), std::forward<A>(args)...);
        std::shared_ptr<std::promise<R>> task_promise = std::make_shared<std::promise<R>>();
        push_task(
            [task_function, task_promise] {
                try {
                    if constexpr(std::is_void_v<R>) {
                        std::invoke(task_function);
                        task_promise->set_value();
                    }
                    else { task_promise->set_value(std::invoke(task_function)); }
                }
                catch(...) {
                    try { task_promise->set_exception(std::current_exception()); }
                    catch(...) {}
                }
            });
        return task_promise->get_future();
    }

    void wait_for_tasks() {
        if(waiting) return;
        waiting = true;
        std::unique_lock tasks_lock(tasks_mutex);
        task_done_cv.wait(tasks_lock, [this] { return (tasks_total == 0); });
        waiting = false;
    }
};

#endif
