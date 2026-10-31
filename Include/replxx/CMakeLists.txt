add_library(replxx STATIC
        ConvertUTF.cpp
        conversion.cxx
        escape.cxx
        history.cxx
        prompt.cxx
        replxx.cxx
        replxx_impl.cxx
        terminal.cxx
        util.cxx
        wcwidth.cpp
        windows.cxx
)

target_compile_definitions(replxx PUBLIC REPLXX_STATIC)

find_package(Threads REQUIRED)

target_link_libraries(replxx PUBLIC Threads::Threads)
