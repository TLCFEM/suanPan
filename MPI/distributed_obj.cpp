#include "Distributed.h"

#include <mpl/mpl.hpp>

class Object : public Distributed {
    mat data{-999, -999};

public:
    Object(int tag)
        : Distributed(tag) {}

    auto generate() {
        if(is_local) data.fill(fill::randn);
    }

    auto gather_data() { return gather(data); }

    auto print() { printf("Object %d on process %d data: %+.6f %+.6f.\n", tag, mpl::environment::comm_world().rank(), data(0), data(1)); }
};

int main() {
    arma_rng::set_seed_random();

    std::vector<Object> objects;

    for(auto i = 0; i < 100; ++i) objects.emplace_back(i);

    for(auto& i : objects) i.generate();

    mpl::irequest_pool requests;
    for(auto& i : objects)
        if(auto req = i.gather_data(); req.has_value()) requests.push(std::move(req).value());
    requests.waitall();

    if(0 == mpl::environment::comm_world().rank())
        for(auto& i : objects) i.print();

    return 0;
}