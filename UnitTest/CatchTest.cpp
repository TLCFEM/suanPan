#include "CatchTest.h"
#include <Toolbox/utility.h>
#include "CatchHeader.h"

void catchtest_main(const int argc, char** argv) {
    for(auto I = 1; I < argc; ++I) {
        constexpr auto t_argv = "";
        bool flag = false;
        if(is_equal(argv[I], "-ctest") || is_equal(argv[I], "--catchtest")) flag = true;
        argv[I] = const_cast<char*>(t_argv);
        if(flag) break;
    }

    Catch::Session().run(argc, argv);
}
