//
// Created by marcel on 26.05.21.
//

#ifndef FICTION_EQUIVALENCE_CHECKING_UTILS_HPP
#define FICTION_EQUIVALENCE_CHECKING_UTILS_HPP

#include <catch2/catch_test_macros.hpp>

#include <fiction/algorithms/verification/equivalence_checking.hpp>

template <typename Spec, typename Impl>
void check_eq(const Spec& spec, const Impl& impl)
{

    fiction::equivalence_checking_stats st{};
    fiction::equivalence_checking(spec, impl, &st);

    CHECK(st.eq != fiction::eq_type::NO);
}

#endif  // FICTION_EQUIVALENCE_CHECKING_UTILS_HPP
