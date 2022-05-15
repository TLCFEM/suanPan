#define CATCH_CONFIG_RUNNER
#include "CatchTest.h"
#include <Toolbox/utility.h>
#include "CatchHeader.h"

int catchtest_main(const int argc, char** argv) {

    for(auto I = 1; I < argc; ++I) if(constexpr auto t_argv = ""; is_equal(argv[I], "-ctest") || is_equal(argv[I], "--catchtest")) argv[I] = const_cast<char*>(t_argv);

    return Catch::Session().run(argc, argv);
}
