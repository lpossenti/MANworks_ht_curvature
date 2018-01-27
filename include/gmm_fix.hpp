#ifndef GMM_FIX_HPP_
#define GMM_FIX_HPP_

#include <gmm/gmm_def.h>
#include <gmm/gmm_sub_vector.h>
#include <gmm/gmm_scaled.h>

namespace gmm {
    template <class V, class S, class SUBI>
    auto
    sub_vector(const scaled_vector_const_ref<V, S> &v, const SUBI &si)
    {
        using tab_ref = tab_ref_with_origin<
            typename scaled_vector_const_ref<V, S>::iterator, V>;
        auto sv = sub_vector(tab_ref(v.begin_, v.end_, v.origin), si);

        return scaled_vector_const_ref<decltype(sv), S>(sv, v.r);
    }


    struct owned_implementation {};

    template <>
    struct principal_orientation_type<owned_implementation> {
        using potype = owned_implementation;
    };

    template <class L1, class L2, class L3>
    inline
    void mult_spec(const L1 &m, const L2 &src, L3 &dst, owned_implementation tag)
    {
        m.mult(src, dst);
    }
} // namespace gmm

#endif // GMM_FIX_HPP_
