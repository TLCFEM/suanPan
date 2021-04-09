#define CATCH_CONFIG_RUNNER
#define CATCH_CONFIG_ENABLE_BENCHMARKING
#include "CatchTest.h"
#include <catch/catch.hpp>
#include <Toolbox/utility.h>

int catchtest_main(const int argc, char** argv) {

	for(auto I = 1; I < argc; ++I) if(const auto t_argv = ""; is_equal(argv[I], "-ctest") || is_equal(argv[I], "--catchtest")) argv[I] = const_cast<char*>(t_argv);

	return Catch::Session().run(argc, argv);
}
