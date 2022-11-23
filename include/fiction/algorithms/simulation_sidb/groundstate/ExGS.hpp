//
// Created by Jan Drewniok on 23.11.22.
//

#ifndef FICTION_EXGS_HPP
#define FICTION_EXGS_HPP

namespace fiction
{
struct physParams
{
    static constexpr float epsilon        = 8.854f * 1E-12f;
    static constexpr float epsilon_screen = 5.6f;
    static constexpr float k              = 1.0f / (4.0f * 3.141592653f * epsilon * epsilon_screen);
    static constexpr float e              = 1.602f * 1E-19f;
    static constexpr float mu             = -0.32f;
    static constexpr float mu_p           = mu - 0.59f;
    static constexpr float tf             = 5.0f * 1E-9f;
    static constexpr float POP_STABILITY_ERR = 1E-6f;
};


namespace detail
{

}



}



#endif  // FICTION_EXGS_HPP
