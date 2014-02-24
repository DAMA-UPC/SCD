/*SCD is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

SCD is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef SCD_TIME_H
#define SCD_TIME_H

#include <common/types.h>
#include <stdlib.h>
#include <sys/time.h>

namespace scd {

    /** @brief Gets the current time in miliseconds.
     *  @return The time in miliseconds.**/
    uint64_t StartClock();

    /** @brief Gets the time elapsed from a given moment.
     *  @param[in] The time specifying the moment from whith to compute the time in miliseconds
     *  @return The time elapsed since the specifyied time in miliseconds.**/
    uint64_t StopClock(uint64_t initTime);

}
#endif
