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

#include <common/time.h>

namespace scd {

    uint64_t StartClock() {
        timeval time;
        gettimeofday(&time, NULL);
        uint64_t initTime = (time.tv_sec * 1000) + (time.tv_usec / 1000);
        return initTime;
    }

    uint64_t StopClock(uint64_t initTime) {
        timeval time;
        gettimeofday(&time, NULL);
        uint64_t endTime = (time.tv_sec * 1000) + (time.tv_usec / 1000);
        return endTime - initTime;
    }
}
