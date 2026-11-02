#include "CatchHeader.h"

#include <Toolbox/thread_pool.hpp>
#include <atomic>
#include <chrono>
#include <numeric>

TEST_CASE("Thread Pool Submit Values", "[Utility.ThreadPool]") {
    thread_pool pool(4);
    REQUIRE(4 == pool.get_thread_count());

    auto f1 = pool.submit([](int a, int b) { return a + b; }, 3, 4);
    auto f2 = pool.submit([](int a) { return a * a; }, 7);
    auto f3 = pool.submit([] { return 42; });

    pool.wait_for_tasks();

    REQUIRE(7 == f1.get());
    REQUIRE(49 == f2.get());
    REQUIRE(42 == f3.get());
}

TEST_CASE("Thread Pool Submit On Default Thread Count", "[Utility.ThreadPool]") {
    thread_pool pool; // auto thread count
    REQUIRE(pool.get_thread_count() > 0);

    auto future = pool.submit([](double x) { return x * 2.; }, 2.5);
    pool.wait_for_tasks();
    REQUIRE(Approx(5.) == future.get());
}

TEST_CASE("Thread Pool Vector Accumulation", "[Utility.ThreadPool]") {
    constexpr auto N = 1000;
    std::vector<int> data(N);
    std::iota(data.begin(), data.end(), 0);

    const int expected = (N - 1) * N / 2;

    thread_pool pool(4);

    std::atomic<int> sum{0};

    pool.push_loop(N, [&](const std::size_t start, const std::size_t end) {
        for(std::size_t I = start; I < end; ++I) sum += data[I];
    });
    pool.wait_for_tasks();

    REQUIRE(expected == sum.load());
}

TEST_CASE("Thread Pool Push Loop Bounds", "[Utility.ThreadPool]") {
    thread_pool pool(3);

    // the 3-argument push_loop passes each block's [start, end) range to the loop
    std::atomic<int> counter{0};
    pool.push_loop(2, 7, [&](const std::size_t start, const std::size_t end) {
        for(std::size_t i = start; i < end; ++i) counter += static_cast<int>(i);
    });
    pool.wait_for_tasks();

    REQUIRE(2 + 3 + 4 + 5 + 6 == counter.load());
}

TEST_CASE("Thread Pool Reversed Push Loop Bounds", "[Utility.ThreadPool]") {
    // push_loop swaps reversed iteration bounds
    thread_pool pool(3);

    std::atomic<int> counter{0};
    pool.push_loop(7, 2, [&](const std::size_t start, const std::size_t end) {
        for(std::size_t i = start; i < end; ++i) ++counter;
    });
    pool.wait_for_tasks();

    REQUIRE(5 == counter.load());
}

TEST_CASE("Thread Pool Submit Exception Propagation", "[Utility.ThreadPool]") {
    thread_pool pool(2);

    auto future = pool.submit([]() -> int { throw std::runtime_error("boom"); });
    pool.wait_for_tasks();

    REQUIRE_THROWS_AS(future.get(), std::runtime_error);
}

TEST_CASE("Thread Pool Void Task", "[Utility.ThreadPool]") {
    thread_pool pool(2);

    std::atomic<bool> done{false};
    auto future = pool.submit([&done] { done = true; });
    pool.wait_for_tasks();
    future.get();

    REQUIRE(done.load());
}

TEST_CASE("Thread Pool Many Tasks", "[Utility.ThreadPool]") {
    thread_pool pool(8);

    constexpr auto N = 500;
    std::vector<std::future<int>> futures;
    futures.reserve(N);

    for(auto I = 0; I < N; ++I) futures.emplace_back(pool.submit([I] { return I * I; }));

    pool.wait_for_tasks();

    for(auto I = 0; I < N; ++I) {
        const auto v = futures[I].get();
        REQUIRE(I * I == v);
    }
}
