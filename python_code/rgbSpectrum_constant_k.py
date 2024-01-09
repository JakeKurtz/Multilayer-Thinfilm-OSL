# ------------------------------------------------------------------------- #
#
#    Copyright (C) 2023 Jake Kurtz
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program. If not, see <https://www.gnu.org/licenses/>.
#
# ------------------------------------------------------------------------- #

# title "Rendering biological iridescences with RGB-based renderers"
# url = https://doi.org/10.1145/1122501.1122506
# doi = 10.1145/1122501.1122506

# The solution for equation (5.7) 
# Solves for the constant k.

from spectral import *

k = 0
for i in range(CIE_SAMPLES):
    k += CIE_D65[i]/CIE_Y_integral * CIE_Y[i] * CIE_STEP

print("k = "+str(k))